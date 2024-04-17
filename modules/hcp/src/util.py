"""Module for fetching files from HCP."""
import sys
from pathlib import Path

from NGPIris.hcp import HCPManager


def fetch(
    *,
    credentials: Path,
    local_path: Path,
    remote_key: str,
) -> Path:
    """Fetches a file from HCP."""
    sys.stdout = open("/dev/null", "w", encoding="utf-8")
    sys.stderr = open("/dev/null", "w", encoding="utf-8")

    hcpm = HCPManager(
        credentials_path=credentials,
        bucket="data",  # FIXME: make this configurable
    )

    hcpm.download_file(
        remote_key,
        local_path=str(local_path),
        callback=False,
        force=True,
    )
    return local_path


def callback(sample, f_idx, logger):
    """Callback for fetching files from HCP."""
    def inner(local_path: Path):
        logger.debug(f"Fetched {local_path.name} from hcp")
        sample.files.insert(f_idx, local_path)

    return inner


def error_callback(sample, logger):
    """Error callback for fetching files from HCP."""
    def inner(exception: Exception):
        logger.error(f"Failed to fetch backup for {sample.id} ({exception})")
        sample.files = []

    return inner
