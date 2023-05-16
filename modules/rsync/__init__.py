import multiprocessing as mp
from functools import partial
from pathlib import Path
from logging import LoggerAdapter
from itertools import groupby
from humanfriendly import parse_size

from cellophane import modules, data, cfg, sge


def _sync_callback(
    logger: LoggerAdapter,
    outputs: list[data.Output],
):
    logger.debug(f"Synced {len(outputs)} outputs")
    for o in outputs:
        dest = o.dest_dir / o.src[0].name
        if dest.exists():
            logger.debug(f"Copied {o.src[0]} to {o.dest_dir}")
        else:
            logger.warning(f"{o.src[0]} is missing from {o.dest_dir}")


def _sync_error_callback(
    code: int,
    logger: LoggerAdapter,
    outputs: list[data.Output],
):
    if len(outputs) == 1:
        logger.error(f"Sync failed for {outputs[0].src} ({code=})")
    else:
        logger.error(f"Sync failed for {len(outputs)} outputs ({code=})")


def _group_by_dest_dir(outputs: list[data.Output]):
    return [
        data.Output(src=[s for o in g for s in o.src], dest_dir=k)
        for k, g in groupby(sorted(outputs), lambda o: o.dest_dir)
    ]


@modules.post_hook(label="Sync Output", condition="complete")
def rsync_results(
    samples: data.Samples,
    logger: LoggerAdapter,
    config: cfg.Config,
    **_,
) -> None:
    if config.rsync.skip:
        logger.info("Skipping output sync")
        return
    elif not any(s.output for s in samples):
        logger.warning("No output to sync")
        return
    else:
        logger.info(f"Syncing output to {config.rsync.base}")

    _outputs = [o for s in samples for o in s.output or []]

    # Split outputs into large files, small files, and directories
    _large_files: list[data.Output] = []
    _small_files: list[data.Output] = []
    _directories: list[data.Output] = []
    for output in _outputs:
        for src in output.src:
            _output = data.Output(src=src, dest_dir=output.dest_dir)
            if src.is_dir():
                _directories.append(_output)
            elif src.stat().st_size > parse_size(config.rsync.large_file_threshold):
                _large_files.append(_output)
            else:
                _small_files.append(_output)

    # Merge outputs with the same destination directory
    _large_files = _group_by_dest_dir(_large_files)
    _small_files = _group_by_dest_dir(_small_files)
    _directories = _group_by_dest_dir(_directories)

    _procs: list[mp.Process] = []
    for label, category in (
        ("large files", _large_files),
        ("small files", _small_files),
        ("directories", _directories),
    ):
        if category:
            logger.debug(f"Syncing {len(category)} {label}")

        _proc = sge.submit(
            str(Path(__file__).parent / "scripts" / "rsync.sh"),
            queue=config.rsync.sge_queue,
            pe=config.rsync.sge_pe,
            slots=config.rsync.sge_slots,
            name="rsync",
            check=False,
            env={
                "SRC": " ".join(",".join(str(s) for s in o.src) for o in category),
                "DST": " ".join(str(o.dest_dir) for o in category),
            },
            callback=partial(
                _sync_callback,
                logger=logger,
                outputs=category,
            ),
            error_callback=partial(
                _sync_error_callback,
                logger=logger,
                outputs=category,
            ),
        )
        _procs.append(_proc)

    for proc in _procs:
        proc.join()

    logger.info("Finished syncing output")
