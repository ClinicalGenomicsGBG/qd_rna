from copy import deepcopy
from logging import LoggerAdapter
from pathlib import Path

from cellophane import Config, Output, Sample, Samples, pre_hook
from git import InvalidGitRepositoryError, NoSuchPathError, Repo

from modules.slims import SlimsSamples


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


@pre_hook(label="Find Linked", after=["slims_fetch"], before=["hcp_fetch", "slims_derive"])
def get_linked_samples(
    samples: SlimsSamples,
    logger: LoggerAdapter,
    config: Config,
    **_,
) -> Samples:
    logger.debug("Fetching samples from earlier runs")
    criteria = "{base_criteria} and ({link_criteria})".format(
        base_criteria=config.slims.find_criteria,
        link_criteria=" or ".join(
            f"(cntn_id equals {group} "
            f"and cntn_cstm_runTag not_one_of {' '.join(s.meta.run for s in samples)})"
            for group, samples in samples.split(by="id")
        ),
    )
    linked_samples = samples.__class__.from_criteria(criteria=criteria, config=config)
    logger.info(f"Found {len(linked_samples)} linked records")
    return linked_samples | samples


@pre_hook(label="Sample ID", after=["unpack"], before=["start_mail"])
def set_sample_id(
    samples: Samples[Sample],
    logger: LoggerAdapter,
    config: Config,
    workdir: Path,
    **_,
) -> Samples:
    logger.debug("Adding runtag to sample IDs")
    _samples = deepcopy(samples)
    known_dups: set[str] = set()
    for sample in samples:
        dups = [s for s in _samples if s.id == sample.id]
        runs = set(sorted(d.meta.run for d in dups if "run" in d.meta))

        if (n := len(dups)) > 1 and config.merge:
            if sample.id not in known_dups:
                logger.info(f"{n} samples with id {sample.id} will be merged")
            known_dups |= {sample.id}
            sample.id = f"{sample.id}_{sorted(runs)[-1]}" if runs else sample.id
            workdir.mkdir(exist_ok=True)
            merge_file = workdir / f"{sample.id}.merged_runs.txt"
            merge_file.write_text("\n".join(runs))
            samples.output.add(
                Output(src=merge_file, dst=Path(sample.id) / merge_file.name)
            )
        elif n > 1 and not config.merge and "run" in sample.meta:
            sample.id = f"{sample.id}_{sample.meta.run}"
        elif n > 1 and not config.merge:
            if sample.id not in known_dups:
                logger.warning(f"Will merge {n} samples with shared id {sample.id}")
            known_dups |= {sample.id}

        else:
            sample.id = (
                f"{sample.id}_{sample.meta.run}" if "run" in sample.meta else sample.id
            )

    return samples
