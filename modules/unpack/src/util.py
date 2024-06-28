"""Module for fetching files from HCP."""

from copy import copy
from logging import LoggerAdapter
from pathlib import Path
from time import sleep

from cellophane import Cleaner, Sample

from .extractors import Extractor


def callback(
    result: None,
    /,
    extractor: Extractor,
    timeout: int,
    sample: Sample,
    idx: int,
    logger: LoggerAdapter,
    path: Path,
    cleaner: Cleaner,
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
        logger.error(
            f"Extracted files for {path.name} "
            f"not found after {timeout} seconds"
        )

    for extracted_path in extracted_paths:
        logger.debug(f"Extracted {extracted_path.name}")
        sample.files.insert(idx, extracted_path)
        cleaner.register(extracted_path.resolve(), ignore_outside_root=True)

    if path in sample.files:
        sample.files.remove(path)


def error_callback(
    exception: Exception,
    /,
    sample: Sample,
    logger: LoggerAdapter,
    path: Path,
    extractor: Extractor,
    cleaner: Cleaner,
) -> None:
    logger.error(f"Failed to extract {path.name}: {exception!r}")
    sample.fail(f"Failed to extract {path.name}")
    if path in sample.files:
        sample.files.remove(path)
    for extracted_path in extractor.extracted_paths(path):
        cleaner.register(extracted_path.resolve(), ignore_outside_root=True)
