import multiprocessing as mp
from functools import partial
from pathlib import Path
from logging import LoggerAdapter, WARNING, DEBUG
import sys
import os
from humanfriendly import parse_size

from cellophane import modules, data, cfg, sge


def _sync_callback(
    logger: LoggerAdapter,
    outputs: list[data.Output],
):
    logger.debug(f"Synced {len(outputs)} outputs")
    for o in outputs:
        dest = (o.dest_dir / o.src[0].name)
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


@modules.post_hook(label="Output", condition="complete")
def rsync_results(
    samples: data.Samples,
    logger: LoggerAdapter,
    config: cfg.Config,
    **_,
) -> None:

    if not any(s.output for s in samples):
        return

    logger.info(f"Copying output to {config.rsync.base}")

    # Generate a set of unique outputs
    unique_outputs: set[data.Output] = set()
    for sid, output in ((s.id, o) for s in samples for o in s.output or []):

        _dst = output.dest_dir
        
        if not output.dest_dir.is_absolute():
            _dst = (config.rsync.base / _dst).absolute()

        if not _dst.is_relative_to(config.rsync.base):
            logger.error(f"Output is not relative to {config.rsync.base}")
            return
        
        if not output.src:
            logger.warning(f"No source paths for {output.dest_dir}")
        elif [*_dst.glob("*")] and not config.rsync.overwrite:
            logger.error(f"Output directory {output.dest_dir} not empty")
            return
        else:
            
            for src in output.src:
                if not src.exists():
                    logger.warning(f"Source path {src} does not exist")
                else:
                    unique_outputs |= {data.Output(src=src, dest_dir=_dst)}
    
    # Split outputs into chunks
    outputs: list[list[data.Output]] = []
    chunk: list[data.Output] = []
    idx = len(unique_outputs)
    for output in [*unique_outputs]:
        if [*output.src][0].is_dir():
            logger.debug(f"Copying directory {output.src} to {output.dest_dir}")
            _output = data.Output(src=output.src, dest_dir=output.dest_dir.parent)
            outputs.append([output])
            idx -= 1
        elif [*output.src][0].stat().st_size > parse_size(
            config.rsync.large_file_threshold
        ):
            logger.debug(f"Copying large file {output.src} to {output.dest_dir}")
            outputs.append([output])
            idx -= 1
        elif len(chunk) >= config.rsync.files_per_job:
            outputs.append(chunk)
            idx -= len(chunk)
            chunk = []
        else:
            chunk.append(output)

        if len(chunk) == idx:
            logger.debug(f"Copying {len(chunk)} files to {output.dest_dir}")
            outputs.append(chunk)
            break

    for chunk in outputs:
        src, dst = [], []
        for output in chunk:
            src += [output.src[0].absolute()]
            if Path(output.src[0]).is_dir():
                output.dest_dir.parent.mkdir(exist_ok=True)
                dst += [output.dest_dir.absolute()]
            else:
                output.dest_dir.mkdir(parents=True, exist_ok=True)
                dst += [output.dest_dir.absolute() / output.src[0].name]

        logger.debug(f"Syncing chunk: {chunk}")
        sge.submit(
            str(Path(__file__).parent / "scripts" / "rsync.sh"),
            queue=config.rsync.sge_queue,
            pe=config.rsync.sge_pe,
            slots=config.rsync.sge_slots,
            name="rsync",
            env={
                "SRC": " ".join(str(s) for s in src),
                "DST": " ".join(str(d) for d in dst),
            },
            check=True,
            stdout=Path("/home/xdemer/test.out"),
            stderr=Path("/home/xdemer/test.err"),
            callback = partial(
                _sync_callback,
                logger=logger,
                outputs=chunk,
            ),
            error_callback = partial(
                _sync_error_callback,
                logger=logger,
                outputs=chunk,
            ),
            temporary_local_flag_please_delete=True
        )

    logger.debug("Finished syncing output")
