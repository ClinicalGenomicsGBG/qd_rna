from dataclasses import dataclass
from pathlib import Path
from shutil import copy
from cellophane import modules, data, cfg
from logging import LoggerAdapter
from typing import Iterable
from concurrent.futures import ProcessPoolExecutor

@dataclass
class Output:
    """Output dataclass for QD-RNA samples."""
    src: Iterable[Path]
    dest: Path

    def __post_init__(self):
        self.dest = Path(self.dest)
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

@modules.pre_hook(label="Set Sample ID")
def set_sample_id(samples: QDRNASamples, logger: LoggerAdapter, config: cfg.Config, **_):
    logger.debug("Adding Run ID to sample IDs")
    for idx, s in enumerate(samples):
        samples[idx].id = f"{s.id}_{s.run}" if "run" in s and s.run else s.id

@modules.post_hook(label="Copy results")
def copy_results(samples: QDRNASamples, logger: LoggerAdapter, config: cfg.Config, **_):
    logger.info(f"Copying output to {config.results.base}")
    with ProcessPoolExecutor(config.results.parallel) as executor:
        for output in [o for s in samples if s.complete for o in s.output]:
            for src in output.src:
                try:
                    dest = config.results.base / output.dest
                    dest.mkdir(parents=True, exist_ok=config.results.overwrite)
                    logger.debug(f"Copying {src} to {dest}")
                    executor.submit(copy, src, dest)
                except Exception as e:
                    logger.warning(f"Failed to copy {src} to {dest}: {e}")

