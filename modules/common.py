"""Common hooks and utilities for QD-RNA."""

from logging import LoggerAdapter
from pathlib import Path

from attrs import define, field
from cellophane import Config, Sample, Samples, pre_hook
from git import InvalidGitRepositoryError, NoSuchPathError, Repo

from modules.slims import SlimsSamples


@define(slots=False, init=False)
class QDRRNASample(Sample):
    """Sample with run information."""
    run: str | None = field(default="UNSPECIFIED")
    last_run: str | None = field(default="UNSPECIFIED")

    @staticmethod
    @Sample.merge.register("run")
    @Sample.merge.register("last_run")
    def _merge(this, that):
        return this or that


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
    before=["hcp_fetch", "slims_derive", "set_sample_id"],
)
def get_linked_samples(
    samples: SlimsSamples,
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
        base_criteria=config.slims.find_criteria,
        link_criteria=" or ".join(
            f"(cntn_id equals {group} "
            f"and cntn_cstm_runTag not_one_of {' '.join(s.run for s in samples)})"
            for group, samples in samples.split(by="id")
        ),
    )
    linked_samples = samples.__class__.from_criteria(criteria=criteria, config=config)
    logger.info(f"Found {len(linked_samples)} linked records")
    return linked_samples | samples


@pre_hook(label="Update Most Recent Run", before=["start_mail"])
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
    logger.debug("Updating most recent run")
    for _, group in samples.split(by="id"):
        # FIXME: Ideally we should sort by the run date.
        # This information is currently  not available. A possible solution
        # would be to allow mapping from parent records in the SLIMS module,
        # but this requires a major refactor of the module.
        latest = max(group, key=lambda s: s.run).run
        for sample in group:
            sample.last_run = latest

    return samples
