"""Utility functions for rsync module."""

from copy import copy
from logging import LoggerAdapter
from pathlib import Path
from time import sleep


def sync_callback(
    result: None,
    /,
    logger: LoggerAdapter,
    manifest: list[tuple[str, str]],
    timeout: int,
):
    """Callback function for rsync_results. Waits for files to become available."""
    del result  # Unused
    for src, dst in manifest:
        if not Path(dst).exists():
            logger.debug(f"Waiting {timeout} seconds for {dst} to become available")
        _timeout = copy(timeout)
        while not (available := Path(dst).exists()) and (_timeout := _timeout - 1) > 0:
            sleep(1)
        if available:
            logger.debug(f"Copied {src} -> {dst}")
        else:
            logger.warning(f"{dst} is missing")
