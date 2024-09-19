"""Module for fetching files from HCP."""

from functools import partial
from logging import LoggerAdapter
from pathlib import Path

from cellophane import Cleaner, Config, Executor, Samples, pre_hook
from mpire.async_result import AsyncResult

from .extractors import Extractor, PetageneExtractor, SpringExtractor
from .util import callback, error_callback

EXTRACORS: dict[str, Extractor] = {
    ".fasterq": PetageneExtractor(),
    ".spring": SpringExtractor(),
}


@pre_hook(label="unpack", after=["hcp_fetch"])
def unpack(
    samples: Samples,
    config: Config,
    logger: LoggerAdapter,
    executor: Executor,
    cleaner: Cleaner,
    workdir: Path,
    **_,
) -> Samples:
    """Extract petagene fasterq files."""
    results: list[AsyncResult] = []
    for sample, idx, path, extractor in (
        (s, i, p, EXTRACORS[Path(p).suffix])
        for s in samples
        for i, p in enumerate(s.files)
        if Path(p).suffix in EXTRACORS
    ):
        (workdir / "unpack").mkdir(parents=True, exist_ok=True)
        if result := extractor.extract(
            logger=logger,
            compressed_path=path,
            config=config,
            executor=executor,
            workdir=workdir / "unpack",
            callback=partial(
                callback,
                extractor=extractor,
                timeout=config.unpack.timeout,
                sample=sample,
                idx=idx,
                logger=logger,
                path=path,
                cleaner=cleaner,
                workdir=workdir / "unpack",
            ),
            error_callback=partial(
                error_callback,
                sample=sample,
                logger=logger,
                path=path,
                extractor=extractor,
                cleaner=cleaner,
                workdir=workdir / "unpack",
            ),
        ):
            results.append(result)

    executor.wait()
    return samples
