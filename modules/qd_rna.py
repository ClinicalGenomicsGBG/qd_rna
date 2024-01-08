from pathlib import Path
from cellophane import cfg, data, modules
from logging import LoggerAdapter
from copy import deepcopy


@modules.pre_hook(label="Sample ID", after=["unpack"], before=["start_mail"])
def set_sample_id(
    samples: data.Samples[data.Sample],
    logger: LoggerAdapter,
    config: cfg.Config,
    workdir: Path,
    **_,
) -> data.Samples:
    logger.debug("Adding Run ID to sample IDs")
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
            merge_file = workdir / f"{sample.id}.merged_runs.txt"
            merge_file.write_text("\n".join(runs))
            sample.output += [data.Output(src=merge_file, dst=Path(sample.id).merged_runs.txt)]
        elif n > 1 and not config.merge and "run" in sample.meta:
            sample.id = f"{sample.id}_{sample.meta.run}"
        elif n > 1 and not config.merge and "run" not in sample.meta:
            if sample.id not in known_dups:
                logger.warning(f"Will merge {n} samples with shared id {sample.id}")
            known_dups |= {sample.id}

        else:
            sample.id = f"{sample.id}_{sample.meta.run}" if "run" in sample.meta else sample.id

    return samples
