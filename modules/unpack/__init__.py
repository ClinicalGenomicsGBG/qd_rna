"""Module for fetching files from HCP."""

from copy import copy
from functools import partial
from logging import LoggerAdapter
from pathlib import Path
from time import sleep

from cellophane import cfg, data, executors, modules
from mpire.async_result import AsyncResult

from .src.extractors import Extractor, PetageneExtractor, SpringExtractor

extractors: dict[str, Extractor] = {
    ".fasterq": PetageneExtractor(),
    ".spring": SpringExtractor(),
}


def _callback(
    result: None,
    /,
    extractor: Extractor,
    timeout: int,
    sample: data.Sample,
    idx: int,
    logger: LoggerAdapter,
    path: Path,
) -> None:
    del result  # Unused

    if not [*extractor.extracted_paths(path)]:
        logger.debug(f"Waiting up to {timeout} seconds for files to become available")
    _timeout = copy(timeout)
    while (
        not (extracted_paths := [*extractor.extracted_paths(path)])
        and (_timeout := _timeout - 1) > 0
    ):
        sleep(1)

    if not extracted_paths:
        logger.error(f"Failed to extract {path.name}")

    for extracted_path in extracted_paths:
        logger.debug(f"Extracted {extracted_path.name}")
        sample.files.insert(idx, extracted_path)

    if path in sample.files:
        sample.files.remove(path)


def _error_callback(
    exception: Exception,
    /,
    sample: data.Sample,
    logger: LoggerAdapter,
    path: Path,
) -> None:
    logger.error(f"Failed to extract {path.name}: {exception}")
    if path in sample.files:
        sample.files.remove(path)


@modules.pre_hook(label="unpack", after=["hcp_fetch"])
def unpack(
    samples: data.Samples,
    config: cfg.Config,
    logger: LoggerAdapter,
    executor: executors.Executor,
    **_,
) -> data.Samples:
    """Extract petagene fasterq files."""
    results: list[AsyncResult] = []
    for sample, idx, path, extractor in (
        (s, i, p, extractors[Path(p).suffix])
        for s in samples
        for i, p in enumerate(s.files)
        if Path(p).suffix in extractors
    ):
        if result := extractor.extract(
            logger=logger,
            compressed_path=path,
            config=config,
            executor=executor,
            callback=partial(
                _callback,
                extractor=extractor,
                timeout=config.unpack.timeout,
                sample=sample,
                idx=idx,
                logger=logger,
                path=path,
            ),
            error_callback=partial(
                _error_callback,
                sample=sample,
                logger=logger,
                path=path,
            ),
        ):
            results.append(result)

    executor.wait()
    return samples
