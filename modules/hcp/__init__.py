"""Module for fetching files from HCP."""
import sys
from logging import LoggerAdapter
from pathlib import Path

from attrs import Attribute, define, field
from attrs.setters import convert
from cellophane import cfg, data, modules
from mpire import WorkerPool
from NGPIris.hcp import HCPManager


def _fetch(
    *,
    credentials: Path,
    local_path: Path,
    remote_key: str,
) -> str:
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


def _callback(sample, f_idx, logger):
    def inner(local_path: Path):
        logger.debug(f"Fetched {local_path.name} from hcp")
        sample.files.insert(f_idx, local_path)

    return inner


def _error_callback(sample, logger):
    def inner(exception: Exception):
        logger.error(f"Failed to fetch backup for {sample.id} ({exception})")
        sample.files = []

    return inner


@define(slots=False)
class HCPSample(data.Sample):
    """Sample with HCP backup."""

    hcp_remote_keys: set[str] | None = field(
        default=None,
        kw_only=True,
        converter=lambda value: None if value is None else set(value),
        on_setattr=convert,
    )

    @hcp_remote_keys.validator
    def _validate_hcp_remote_keys(
        self,
        attribute: Attribute,
        value: set[str] | None,
    ) -> None:
        if value is None:
            return
        elif any(not isinstance(v, str) for v in value):
            raise TypeError(f"Invalid type for value in {attribute.name}: {value}")
        elif len(value) != len(self.files):
            raise ValueError(
                f"Length mismatch between {attribute.name} and files: "
                f"{len(value)} != {len(self.files)}"
            )


@data.Sample.merge.register("hcp_remote_keys")
def _(this, that) -> set[str] | None:
    return (this or set()) | (that or set()) if not this or that is None else None


@modules.pre_hook(label="HCP", after=["slims_fetch"])
def hcp_fetch(
    samples: data.Samples,
    config: cfg.Config,
    logger: LoggerAdapter,
    **_,
) -> data.Samples:
    """Fetch files from HCP."""
    for sample in samples.with_files:
        logger.debug(f"Found all files for {sample.id} locally")

    if any(k not in config.get("hcp", {}) for k in ["credentials", "fastq_temp"]):
        logger.warning("HCP not configured")
        return samples

    if samples.without_files:
        logger.info(f"fetching {len(samples.without_files)} samples from HCP")

    with WorkerPool(
        n_jobs=config.hcp.parallel,
        use_dill=True,
    ) as pool:
        for sample in samples.without_files:
            sample.files = []
            if sample.hcp_remote_keys is None:
                logger.warning(f"No backup for {sample.id}")
                continue

            for f_idx, remote_key in enumerate(sample.hcp_remote_keys):
                local_path = config.hcp.fastq_temp / Path(remote_key).name

                if local_path.exists():
                    sample.files.insert(f_idx, local_path)
                    logger.debug(f"Found {local_path.name} locally")
                    continue

                pool.apply_async(
                    _fetch,
                    kwargs={
                        "credentials": config.hcp.credentials,
                        "local_path": local_path,
                        "remote_key": remote_key,
                    },
                    callback=_callback(sample, f_idx, logger),
                    error_callback=_error_callback(sample, logger),
                )

        pool.join()

    return samples
