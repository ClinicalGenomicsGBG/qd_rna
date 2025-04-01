"""Hooks for syncing output to a local or remote location."""

from functools import partial
from logging import LoggerAdapter
from pathlib import Path

from cellophane import cfg, data, executors, modules
from humanfriendly import parse_size

from .util import sync_callback

ROOT = Path(__file__).parent.parent


@modules.post_hook(label="RSync Output", condition="complete")
def rsync_results(
    samples: data.Samples,
    logger: LoggerAdapter,
    config: cfg.Config,
    workdir: Path,
    executor: executors.Executor,
    **_,
) -> None:
    """Sync output to a remote server."""
    if not samples.output:
        logger.warning("No output to sync")
        return

    logger.info(f"Syncing output to {config.resultdir}")

    _workdir = workdir / "rsync"
    _workdir.mkdir(parents=True, exist_ok=True)

    # Split outputs into large files, small files, and directories
    labels: dict[str, str] = {
        "large": f"large files (>{config.rsync.large_file_threshold})",
        "small": f"small files (<{config.rsync.large_file_threshold})",
        "dir": "directories",
    }
    manifests: dict[str, list[tuple[Path, Path]]] = {
        "large": [],
        "small": [],
        "dir": [],
    }

    for output in samples.output:
        if not output.src.exists():
            logger.warning(f"{output.src} does not exist")
            continue
        elif output.dst.exists() and not config.rsync.overwrite:
            logger.warning(f"{output.dst} already exists")
            continue
        elif not output.dst.is_relative_to(config.resultdir):
            logger.warning(f"{output.dst} is outside {config.resultdir}")
            continue
        output.dst.parent.mkdir(parents=True, exist_ok=True)
        if output.src.is_dir():
            manifests["dir"] += [(output.src, output.dst)]
        elif output.src.stat().st_size > parse_size(config.rsync.large_file_threshold):
            manifests["large"] += [(output.src, output.dst)]
        else:
            manifests["small"] += [(output.src, output.dst)]

    for type_, manifest in manifests.items():
        if manifest:
            logger.info(f"Syncing {len(manifest)} {labels[type_]}")
            manifest_path = _workdir / f"rsync.{type_}.manifest"
            with open(manifest_path, "w", encoding="utf-8") as m:
                m.writelines(
                    [
                        f"{src.absolute()}{'/' if src.is_dir() else ''} "
                        f"{dst.absolute()}\n"
                        for src, dst in manifest
                    ]
                )
            executor.submit(
                str(ROOT / "scripts" / "rsync.sh"),
                name="rsync",
                env={"MANIFEST": str(manifest_path.absolute())},
                workdir=_workdir,
                callback=partial(
                    sync_callback,
                    logger=logger,
                    manifest=manifest,
                    timeout=config.rsync.timeout,
                ),
            )
    executor.wait()
    logger.info("Finished syncing output")
    return samples
