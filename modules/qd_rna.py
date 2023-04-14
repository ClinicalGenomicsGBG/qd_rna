from dataclasses import dataclass
from pathlib import Path
from shutil import copy, copytree
from cellophane import modules, data, cfg
from logging import LoggerAdapter
from typing import Iterable, Optional
from copy import deepcopy


@dataclass
class Output:
    """Output dataclass for QD-RNA samples."""

    src: list[Path]
    dest_dir: Path
    dest_name: Optional[str] = None

    def __post_init__(self):
        self.dest_dir = Path(self.dest_dir)
        if not isinstance(self.src, Iterable):
            self.src = [Path(self.src)]
        else:
            self.src = [Path(s) for s in self.src]


class QDRNASample:
    """QD-RNA sample class."""

    output: list[Output] = []


class QDRNASamples(data.Mixin, sample_mixin=QDRNASample):
    """Mixin for QD-RNA samples."""

    pass


@modules.pre_hook(label="Sample ID", after="all")
def set_sample_id(
    samples: data.Samples[data.Sample],
    logger: LoggerAdapter,
    config: cfg.Config,
    **_,
) -> data.Samples:
    logger.debug("Adding Run ID to sample IDs")
    _samples = deepcopy(samples)
    known_dups: set[str] = set()
    for sample in samples:
        dups = [s for s in _samples if s.id == sample.id]
        runs = set(d.run for d in dups  if "run" in d)
        runs = sorted(runs)

        if (n := len(dups)) > 1 and config.merge:
            if sample.id not in known_dups:
                logger.info(f"{n} samples with id {sample.id} will be merged")
            known_dups |= {sample.id}
            sample.id = f"{sample.id}_{sorted(runs)[-1]}" if runs else sample.id
            merge_file = config.outdir / f"{sample.id}.merged_runs.txt"
            merge_file.write_text("\n".join(runs))
            sample.output += [Output(src=[merge_file], dest_dir=Path(sample.id))]
        elif n > 1 and not config.merge and "run" in sample:
            sample.id = f"{sample.id}_{sample.run}"
        elif n > 1 and not config.merge and "run" not in sample:
            if sample.id not in known_dups:
                logger.warning(f"Will merge {n} samples with shared id {sample.id}")
            known_dups |= {sample.id}

        else:
            sample.id = f"{sample.id}_{sample.run}" if "run" in sample else sample.id

    return samples


@modules.post_hook(label="Results", condition="complete")
def copy_results(
    samples: data.Samples,
    logger: LoggerAdapter,
    config: cfg.Config,
    **_,
) -> None:
    logger.info(f"Copying output to {config.results.base}")
    outputs: list[Output] = [o for s in samples for o in s.output]

    for dest_dir in set(config.results.base / o.dest_dir for o in outputs):
        if [*dest_dir.glob("*")] and not config.results.overwrite:
            logger.error(f"Output directory {dest_dir} not empty")
            return

    for output in outputs:
        if not output.src:
            logger.warning(f"No source paths for {output.dest_dir}")
        else:
            for src in [o for o in output.src if not o.exists()]:
                logger.warning(f"Source path {src} does not exist")
                output.src.remove(src)

    for dest_dir, dest_name, src in set([
        (config.results.base / o.dest_dir, o.dest_name, s)
        for o in outputs
        for s in o.src
    ]):
        try:
            if Path(src).is_dir():
                logger.debug(f"Copying directory {src} to {dest_dir}")
                copytree(src, dest_dir, dirs_exist_ok=True)
            else:
                dest_dir.mkdir(parents=True, exist_ok=True)
                dest_path = dest_dir / (dest_name or src.name)
                logger.debug(f"Copying file {src} to {dest_path}")
                copy(src, dest_path)
        except Exception as e:
            logger.warning(f"Failed to copy {src}: {e}")
