from concurrent.futures import ProcessPoolExecutor, Future
from functools import partial
from pathlib import Path
from logging import LoggerAdapter, WARNING, DEBUG
import sys
import os
from humanfriendly import parse_size

from cellophane import modules, data, cfg, sge


def _sync_output(
    *,
    outputs: list[data.Output],
    config: cfg.Config,
) -> None:
    sys.stdout = open(os.devnull, "w", encoding="utf-8")
    sys.stderr = open(os.devnull, "w", encoding="utf-8")

    src, dst = [], []
    for output in outputs:
        src += [output.src[0].absolute()]
        if Path(output.src[0]).is_dir():
            dst += [output.dest_dir.absolute()]
        else:
            output.dest_dir.mkdir(parents=True, exist_ok=True)
            dst += [output.dest_dir.absolute() / output.src[0].name]

    sge.submit(
        str(Path(__file__).parent / "scripts" / "rsync.sh"),
        queue=config.unpack.sge_queue,
        pe=config.unpack.sge_pe,
        slots=config.unpack.sge_slots,
        name="rsync",
        env={
            "SRC": " ".join(str(s) for s in src),
            "DST": " ".join(str(d) for d in dst),
        },
        check=True,
    )


def _sync_callback(
    future: Future,
    /,
    logger: LoggerAdapter,
    outputs: list[data.Output],
):
    exception = future.exception()

    if exception is None:
        logger.debug(f"Synced {len(outputs)} outputs")
    elif len(outputs) == 1:
        logger.error(f"Sync failed for {outputs[0].src} ({exception})")
    else:
        logger.error(f"Sync failed for {len(outputs)} outputs ({exception})")

    for o in outputs:
        dest = (o.dest_dir / o.src[0].name)
        if dest.exists() and dest.stat().st_size == o.src[0].stat().st_size:
            logger.debug(f"Copied {o.src[0]} to {o.dest_dir}")
        else:
            logger.log(
                level=WARNING if exception is None else DEBUG,
                msg=f"{o.src[0]} is missing from {o.dest_dir}",
            )


@modules.post_hook(label="Output", condition="complete")
def rsync_results(
    samples: data.Samples,
    logger: LoggerAdapter,
    config: cfg.Config,
    **_,
) -> None:
    logger.info(f"Copying output to {config.rsync.base}")

    # Generate a set of unique outputs
    unique_outputs: set[data.Output] = set()
    for sid, output in ((s.id, o) for s in samples for o in s.output):
        if not output.src:
            logger.warning(f"No source paths for {output.dest_dir}")
        elif [*output.dest_dir.glob("*")] and not config.rsync.overwrite:
            logger.error(f"Output directory {output.dest_dir} not empty")
            return
        else:
            for src in output.src:
                if not src.exists():
                    logger.warning(f"Source path {src} does not exist")
                else:
                    output.parent_id = sid
                    unique_outputs |= {
                        data.Output(src, config.rsync.base / output.dest_dir)
                    }

    # Split outputs into chunks
    outputs: list[list[data.Output]] = []
    chunk: list[data.Output] = []
    for output in unique_outputs:
        if output.src[0].is_dir():
            logger.debug(f"Copying directory {output.src} to {output.dest_dir}")
            outputs.append([output])
        elif output.src[0].stat().st_size > parse_size(
            config.rsync.large_file_threshold
        ):
            logger.debug(f"Copying large file {output.src} to {output.dest_dir}")
            outputs.append([output])
        elif len(chunk) >= config.rsync.files_per_job:
            logger.debug(f"Copying {len(chunk)} files to {output.dest_dir}")
            outputs.append(chunk)
            chunk = []
        else:
            chunk.append(output)

    with ProcessPoolExecutor() as pool:
        for chunk in outputs:
            pool.submit(
                _sync_output,
                outputs=chunk,
                config=config,
            ).add_done_callback(
                partial(
                    _sync_callback,
                    logger=logger,
                    outputs=chunk,
                )
            )
