"""Module for fetching files from HCP."""

import sys
from concurrent.futures import Future, ProcessPoolExecutor
from functools import partial
from logging import LoggerAdapter
from pathlib import Path
from typing import Sequence

from attrs import define, field
from cellophane import cfg, data, modules
from NGPIris import hcp


@define(slots=False, init=False)
class HCPSample(data.Sample):
    """Sample with HCP backup."""

    backup: list[str] | None = field(default=None)

    @backup.validator
    def validate_backup(self, attribute: str, value: Sequence[str] | None) -> None:
        if not (
            value is None
            or (isinstance(value, Sequence) and all(isinstance(v, str) for v in value))
        ):
            raise ValueError(f"Invalid {attribute} value: {value}")


def _fetch(
    config: cfg.Config,
    local_path: Path,
    remote_key: str | None = None,
) -> None:
    sys.stdout = open(
        config.logdir / f"iris.{local_path.name}.out", "w", encoding="utf-8"
    )
    sys.stderr = open(
        config.logdir / f"iris.{local_path.name}.err", "w", encoding="utf-8"
    )

    hcpm = hcp.HCPManager(
        credentials_path=config.iris.credentials,
        bucket="data",  # FIXME: make this configurable
    )

    if remote_key is None:
        search_path = local_path
        while Path(search_path).suffix:
            search_path = search_path.with_suffix("")

        result = hcpm.search_objects(search_path.name)
        if result is not None and len(result) == 1:
            remote_key = result[0].key
        else:
            raise ValueError(f"Could not find remote key for {local_path}")
    if not local_path.exists():
        hcpm.download_file(
            remote_key,
            local_path=str(local_path),
            callback=False,
            force=True,
        )


def _fetch_callback(
    future: Future,
    /,
    logger: LoggerAdapter,
    samples: data.Samples,
    local_path: Path,
    s_idx: int,
    f_idx: int,
):
    if (exception := future.exception()) is not None:
        logger.error(f"Failed to fetch {local_path} from HCP ({exception})")
        samples[s_idx].files[f_idx] = None
    else:
        logger.debug(f"Fetched {local_path}")
        samples[s_idx].files[f_idx] = str(local_path)


@modules.pre_hook(label="HCP", after=["slims_fetch"])
def hcp_fetch(
    samples: data.Samples,
    config: cfg.Config,
    logger: LoggerAdapter,
    **_,
) -> data.Samples:
    """Fetch files from HCP."""

    # FIXME: add 'skip' option
    with ProcessPoolExecutor(config.iris.parallel) as pool:
        for s_idx, sample in enumerate(samples):
            if all(
                file is not None and Path(file).exists() for file in sample.files or []
            ):
                logger.debug(f"Files found for {sample.id} {sample.files}")

            elif sample.backup:
                logger.info(f"Fetching files for {sample.id} from HCP")
                if "backup" not in sample:
                    logger.warning(f"Remote key not found for {sample.id}, will search")
                    sample.backup = [None] * len(sample.files)

                for f_idx, local_key in enumerate(sample.files):
                    remote_key: str | None = sample.backup[f_idx]
                    _local_key = local_key or remote_key or f"{sample.id}_{f_idx}"
                    local_path = config.iris.fastq_temp / Path(_local_key).name

                    pool.submit(
                        _fetch,
                        config=config,
                        local_path=local_path,
                        remote_key=remote_key,
                    ).add_done_callback(
                        partial(
                            _fetch_callback,
                            logger=logger,
                            samples=samples,
                            local_path=local_path,
                            s_idx=s_idx,
                            f_idx=f_idx,
                        )
                    )
            else:
                logger.warning(f"Unable to fetch files for {sample.id} from HCP")

    return samples.__class__([s for s in samples if s is not None])
