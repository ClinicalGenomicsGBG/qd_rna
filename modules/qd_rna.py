from dataclasses import dataclass
from pathlib import Path
from shutil import copy
from cellophane import modules, data, cfg
from logging import LoggerAdapter
from typing import Iterable, Optional
from concurrent.futures import ProcessPoolExecutor


@dataclass
class Output:
    """Output dataclass for QD-RNA samples."""

    src: Iterable[Path]
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


@modules.pre_hook(label="Sample ID")
def set_sample_id(
    samples: data.Samples,
    logger: LoggerAdapter,
    **_,
):
    logger.debug("Adding Run ID to sample IDs")
    for idx, s in enumerate(samples):
        samples[idx].id = f"{s.id}_{s.run}" if "run" in s and s.run else s.id


@modules.post_hook(label="Results")
def copy_results(
    samples: data.Samples,
    logger: LoggerAdapter,
    config: cfg.Config,
    **_,
):
    logger.info(f"Copying output to {config.results.base}")
    with ProcessPoolExecutor(config.results.parallel) as executor:
        outputs = [o for s in samples for o in s.output]

        for dest_dir in set(config.results.base / o.dest_dir for o in outputs):
            if dest_dir.glob("*") and not config.results.overwrite:
                logger.error(f"Output directory {dest_dir} not empty")
                return
            else:
                dest_dir.mkdir(parents=True, exist_ok=True)

        for dest_dir, dest_name, src in [
            (o.dest_dir, o.dest_name, s) for o in outputs for s in o.src
        ]:
            try:
                dest = config.results.base / dest_dir / (dest_name or src.name)
                logger.debug(f"Copying {src} to {dest}")
                executor.submit(copy, src, dest)
            except Exception as e:
                logger.warning(f"Failed to copy {src}: {e}")
