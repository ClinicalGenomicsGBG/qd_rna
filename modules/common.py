"""Common hooks and utilities for QD-RNA."""

from functools import partial
from logging import LoggerAdapter
from pathlib import Path
from time import sleep
from typing import Literal
from warnings import warn

from attrs import Attribute, define, field
from attrs.setters import convert, validate
from cellophane import Cleaner, Config, Executor, Sample, Samples, post_hook, pre_hook
from git import InvalidGitRepositoryError, NoSuchPathError, Repo
from humanfriendly.text import pluralize


def _int_or_none(value: str) -> int | None:
    if value is not None:
        try:
            return int(value)
        except Exception:  # pylint: disable=broad-except
            warn(f"Failed to convert {value} to int")


@define(slots=False, init=False)
class QDRRNASample(Sample):
    """Sample with run information."""

    run: str | None = field(default="UNSPECIFIED")
    reads: int | None = field(
        default=None,
        converter=_int_or_none,
        on_setattr=convert,
    )
    last_run: str | None = field(default="UNSPECIFIED")
    slims_state: Literal["running", "complete", "error"] = field(
        default="running",
        on_setattr=validate,
    )

    @staticmethod
    @Sample.merge.register("run")
    @Sample.merge.register("last_run")
    def _merge(this, that):
        return this or that

    @staticmethod
    @Sample.merge.register("reads")
    def _merge_reads(this, that):
        if this != that:
            warn(f"Conflicting read counts: {this} != {that}")
            return None
        return this

    @staticmethod
    @Sample.merge.register("slims_state")
    def _merge_state(this, that):
        return (
            "error"
            if "error" in (this, that)
            else "running"
            if "running" in (this, that)
            else "complete"
        )

    @slims_state.validator
    def _validate_state(
        self,
        attribute: Attribute,
        value: Literal["running", "complete", "error"],
    ) -> None:
        del attribute  # Unused
        if value not in ("running", "complete", "error"):
            raise ValueError(f"Invalid state: {value}")


def _fetch_nf_core(
    name: str,
    tag: str,
    url: str,
    path: Path,
    logger: LoggerAdapter,
) -> None:
    try:
        repo = Repo(path)
        if (current := repo.git.tag(points_at="HEAD")) != tag:
            logger.info(f"Updating {name} from {current or repo.head.commit} to {tag}")
            repo.create_head(tag, repo.tags[tag])
    except (NoSuchPathError, InvalidGitRepositoryError):
        logger.info(f"Fetching {name}@{tag}")
        path.mkdir(parents=True, exist_ok=True)
        Repo.clone_from(url, path, branch=tag)
    except Exception as exc:  # pylint: disable=broad-except
        logger.error(f"Failed to fetch {name}@{tag}", exc_info=exc)
        raise SystemExit(1) from exc


@pre_hook(label="Fetch nf-core pipelines", before=["all"])
def fetch_nfcore_pipelines(
    root: Path,
    logger: LoggerAdapter,
    config: Config,
    **_,
) -> None:
    """Fetch nf-core pipelines."""
    _fetch_nf_core(
        name="nf-core/rnaseq",
        path=root / "dependencies" / "nf-core" / "rnaseq",
        tag=config.rnaseq.nf_tag,
        url=config.rnaseq.nf_url,
        logger=logger,
    )
    _fetch_nf_core(
        name="nf-core/rnafusion",
        path=root / "dependencies" / "nf-core" / "rnafusion",
        tag=config.rnafusion.nf_tag,
        url=config.rnafusion.nf_url,
        logger=logger,
    )


def nf_config(template, location, include: Path | None = None, **kwargs):
    """Write nextflow config."""
    with open(location, "w", encoding="utf-8") as f:
        if include is not None:
            f.write(f"includeConfig '{include}'\n\n")
        f.write(template.format(**kwargs))


@pre_hook(
    label="Find Linked",
    after=["slims_fetch"],
    before=["hcp_fetch", "slims_sync_pre", "set_sample_id"],
)
def get_linked_samples(
    samples: Samples,
    logger: LoggerAdapter,
    config: Config,
    **_,
) -> Samples:
    """Find linked samples from earlier runs.

    Samples are linked if they have the same ID but different run tags (i.e. the same
    sample was sequenced multiple times). The find_criteria is used to filter samples
    that should be picked up by QD-RNA.
    """
    logger.debug("Fetching samples from earlier runs")
    criteria = "{base_criteria} and ({link_criteria})".format(
        base_criteria=config.slims.criteria,
        link_criteria=" or ".join(
            f"(cntn_id equals {group} "
            f"and cntn_cstm_runTag not_one_of {' '.join(s.run for s in samples)})"
            for group, samples in samples.split(by="id")
        ),
    )
    linked_samples = samples.__class__.from_criteria(criteria=criteria, config=config)
    logger.info(f"Found {len(linked_samples)} linked records")
    return linked_samples | samples


@post_hook(label="Set SLIMS State", before=["slims_sync_post"])
def set_slims_state(
    samples: Samples,
    logger: LoggerAdapter,
    **_,
) -> Samples:
    """Set SLIMS state based on the presence of a run tag."""
    logger.info("Updating state for SLIMS bioinformatics")
    for sample in samples:
        sample.slims_state = "error" if sample.failed else "complete"


@pre_hook(label="Update Most Recent Run", before=["start_mail"], after="all")
def update_most_recent_run(
    samples: Samples,
    logger: LoggerAdapter,
    **_,
) -> Samples:
    """Update most recent run.

    This is used by QD-RNA to name the output directory of merged samples to
    ID + most recent run. This assumes that sorting the run alphanumerically
    will give the most recent run.
    """
    for id_, group in samples.split(by="id"):
        # FIXME: Ideally we should sort by the run date.
        # This information is currently not available.
        # A possible solution would be to allow mapping from
        # parent records in the SLIMS module, but this requires
        # a major refactor of the module.
        latest = max(group, key=lambda s: s.run).run
        logger.debug(f"Setting most recent run for {id_} to {latest}")
        for sample in group:
            sample.last_run = latest

    return samples


def _subsample_callback(
    result: None,
    /,
    sample: Sample,
    files: tuple[Path, Path],
    logger: LoggerAdapter,
    cleaner: Cleaner,
) -> None:
    del result  # Unused
    if not all(f.exists for f in files):
        logger.debug("Waiting up to 60 seconds for files to become available")
        _timeout = 60
        while not all(f.exists() for f in files) and (_timeout := _timeout - 1) > 0:
            sleep(1)

        if not all(f.exists() for f in files):
            logger.error(f"Subsampled files for {sample.id} not found after 60 seconds")
            sample.fail("Subsampled files not found")
            return

    logger.debug(f"Subsampling finished for {sample.id}")
    sample.files = [*files]
    for f in files:
        cleaner.register(f.resolve())


def _subsample_error_callback(
    exception: Exception,
    /,
    sample: Sample,
    files: tuple[Path, Path],
    logger: LoggerAdapter,
    cleaner: Cleaner,
) -> None:
    logger.error(f"Failed to subsample {sample.id}: {exception}")
    sample.fail(f"Failed to subsample: {exception}")
    for f in files:
        cleaner.register(f.resolve())


@pre_hook(label="Subsample", after=["all", "start_mail"])
def subsample(
    samples: Samples[Sample],
    config: Config,
    logger: LoggerAdapter,
    executor: Executor,
    root: Path,
    workdir: Path,
    cleaner: Cleaner,
    **_,
) -> Samples:
    """Subsample input FASTQs.

    This hook is used to subsample input FASTQs to a fixed number of reads.
    """
    if not config.subsample.target:
        return samples

    for id_, group in samples.split(by="id"):
        if any(sample.reads is None for sample in group):
            logger.info(f"Unknown number of reads for {id_} - not subsampling")
            continue
        elif (total := sum(sample.reads for sample in group)) < config.subsample.target:
            logger.info(f"Too few reads for {id_} - not subsampling")
            continue

        frac = config.subsample.target / total
        n_files = sum(len(sample.files) for sample in group)
        logger.info(
            f"Subsampling {pluralize(n_files, 'file', 'files')} "
            f"for sample {id_} to ~{config.subsample.target} reads "
            f"({frac:.2%})"
        )

        (workdir / "subsample").mkdir(exist_ok=True, parents=True)
        for sample in group:
            subsample_files = (
                workdir / "subsample" / f"{sample.id}_{sample.last_run}_R1.fastq.gz",
                workdir / "subsample" / f"{sample.id}_{sample.last_run}_R2.fastq.gz",
            )
            executor.submit(
                str(root / "scripts" / "common_subsample.sh"),
                name=f"subsample_{id_}",
                workdir=workdir,
                cpus=config.subsample.threads,
                env={
                    "_SUBSAMPLE_INIT": config.qlucore.subsample.init,
                    "_SUBSAMPLE_THREADS": config.qlucore.subsample.threads,
                    "_SUBSAMPLE_FRAC": frac,
                    "_SUBSAMPLE_INPUT_FQ1": sample.files[0],
                    "_SUBSAMPLE_INPUT_FQ2": sample.files[1],
                    "_SUBSAMPLE_OUTPUT_FQ1": subsample_files[0],
                    "_SUBSAMPLE_OUTPUT_FQ2": subsample_files[1],
                },
                callback=partial(
                    _subsample_callback,
                    logger=logger,
                    sample=sample,
                    files=subsample_files,
                    cleaner=cleaner,
                ),
                error_callback=partial(
                    _subsample_error_callback,
                    logger=logger,
                    sample=sample,
                    files=subsample_files,
                    cleaner=cleaner,
                ),
                conda_spec={
                    "dependencies": ["fq >= 0.11.0 < 0.12.0"],
                    "channels": ["bioconda"],
                },
            )
    executor.wait()
    return samples
