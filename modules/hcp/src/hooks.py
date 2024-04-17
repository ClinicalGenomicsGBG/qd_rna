"""Hooks for fetching files from HCP."""

from logging import LoggerAdapter
from pathlib import Path

from cellophane import Config, Samples, pre_hook
from mpire import WorkerPool

from .util import callback, error_callback, fetch


@pre_hook(label="HCP", after=["slims_fetch"])
def hcp_fetch(
    samples: Samples,
    config: Config,
    logger: LoggerAdapter,
    **_,
) -> Samples:
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
                    fetch,
                    kwargs={
                        "credentials": config.hcp.credentials,
                        "local_path": local_path,
                        "remote_key": remote_key,
                    },
                    callback=callback(sample, f_idx, logger),
                    error_callback=error_callback(sample, logger),
                )

        pool.join()

    return samples
