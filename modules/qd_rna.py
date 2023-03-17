from dataclasses import dataclass
from pathlib import Path
from shutil import copyfile
from cellophane import modules, data, cfg
from logging import LoggerAdapter
from typing import Iterable

@dataclass
class Output:
    """Output dataclass for QD-RNA samples."""
    src: Iterable[Path]
    dest: Path

    def __post_init__(self):
        self.destination = Path(self.destination)
        if not isinstance(self.source, Iterable):
            self.source = [Path(self.source)]
        else:
            self.source = [Path(s) for s in self.source]


class QDRNASamples(data.Mixin):
    """Mixin for QD-RNA samples."""
    output: list[Output] = []


@modules.pre_hook(label="QD-RNA Mixin")
def qd_rna_mixin(samples: data.Samples, **_):
    samples.add_mixin(QDRNASamples)
    return samples

@modules.post_hook(label="Copy results")
def copy_webstore(samples: QDRNASamples, logger: LoggerAdapter, config: cfg.Config, **_):
    logger.info("Copying to config.results_base")
    for output in samples.output:
        for src in output.src:
            dest = config.results_base / output.dest
            logger.debug(f"Copying {src} to {dest}")
            copyfile(src, dest)
