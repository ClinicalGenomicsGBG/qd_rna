from copy import deepcopy
from logging import LoggerAdapter
from pathlib import Path

from cellophane import Config, Output, Sample, Samples, pre_hook

from modules.slims import SlimsSamples


@pre_hook(label="Find Linked", after=["slims_fetch"], before=["hcp_fetch"])
def get_linked_samples(
    samples: SlimsSamples,
    logger: LoggerAdapter,
    config: Config,
    **_,
) -> Samples:
    logger.debug("Fetching samples from earlier runs")
    criteria = "{base_criteria} and ({link_criteria})".format(
        base_criteria=config.slims.find_criteria,
        link_criteria="or".join(
            f"(cntn_id equals {sample.id} "
            f"and cntn_cstm_runTag not_equals {sample.meta.run})"
            for sample in samples
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
            samples.output.add(Output(src=merge_file, dst=Path(sample.id) / merge_file.name))
        elif n > 1 and not config.merge and "run" in sample.meta:
            sample.id = f"{sample.id}_{sample.meta.run}"
        elif n > 1 and not config.merge and "run" not in sample.meta:
            if sample.id not in known_dups:
                logger.warning(f"Will merge {n} samples with shared id {sample.id}")
            known_dups |= {sample.id}

        else:
            sample.id = (
                f"{sample.id}_{sample.meta.run}" if "run" in sample.meta else sample.id
            )

    return samples
