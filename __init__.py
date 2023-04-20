from concurrent.futures import ProcessPoolExecutor, Future
from dataclasses import dataclass
from functools import partial
from pathlib import Path
from logging import LoggerAdapter
from typing import Iterable, Optional
import sys
import os

from cellophane import modules, data, cfg, sge

@dataclass
class Output:
    """Output dataclass for samples."""

    src: list[Path]
    dest_dir: Path
    dest_name: Optional[str] = None

    def __post_init__(self):
        self.dest_dir = Path(self.dest_dir)
        if not isinstance(self.src, Iterable):
            self.src = [Path(self.src)]
        else:
            self.src = [Path(s) for s in self.src]


class OutputSample:
    """Mixin for sample with Output files."""

    output: list[Output] = []


class OutputSamples(data.Mixin, sample_mixin=OutputSample):
    """Mixin for samples with Output files."""

    pass


def _sync_file(
    src: Path,
    dst: Path,
    /,
    config: cfg.Config,
) -> None:
    sys.stdout = open(os.devnull, "w", encoding="utf-8")
    sys.stderr = open(os.devnull, "w", encoding="utf-8")

    return sge.submit(
        str(Path(__file__).parent / "scripts" / f"rsync.sh"),
        queue=config.unpack.sge_queue,
        pe=config.unpack.sge_pe,
        slots=config.unpack.sge_slots,
        name="rsync",
        env={
            "SRC": str(src),
            "DST": str(dst),
        },
        cwd=src.parent,
        check=True,
    )


def _sync_callback(
    future: Future,
    /,
    logger: LoggerAdapter,
    src: Path,
    dest: Path,
):
    if (exception := future.exception()) is not None:
        logger.error(f"Failed to copy {src} to {dest} ({exception})")
    elif dest.exists():
        logger.debug(f"Copied {src} to {dest}")
    else:
        logger.error(f"Copy completed, but {dest} does not exist")


@modules.post_hook(label="Output", condition="complete")
def rsync_results(
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

    with ProcessPoolExecutor() as pool:
        for dest_dir, dest_name, src in set([
            (config.results.base / o.dest_dir, o.dest_name, s)
            for o in outputs
            for s in o.src
        ]):
            if Path(src).is_dir():
                dest = dest_dir
                logger.debug(f"Copying directory {src} to {dest_dir}")
            else:
                dest_dir.mkdir(parents=True, exist_ok=True)
                dest = dest_dir / (dest_name or src.name)
                logger.debug(f"Copying file {src} to {dest}")
                
            pool.submit(
                _sync_file,
                src,
                dest,
                config=config,
            ).add_done_callback(
                partial(
                    _sync_callback,
                    logger=logger,
                    src=src,
                    dest=dest,
                )
            )


