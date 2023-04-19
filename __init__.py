"""Module for fetching files from HCP."""

from concurrent.futures import ProcessPoolExecutor, Future
import os
import sys
from functools import partial
from logging import LoggerAdapter
from pathlib import Path
from cellophane import cfg, data, modules, sge


def _extract(
    method: str,
    /,
    compressed_path: Path,
    extract_path: Path,
    config: cfg.Config,
) -> None:
    sys.stdout = open(os.devnull, "w", encoding="utf-8")
    sys.stderr = open(os.devnull, "w", encoding="utf-8")

    match method:
        case "petagene":
            env = dict(
                COMPRESSED_PATH=compressed_path,
                THREADS=config.unpack.sge_slots,
            )
        case "spring":
            env = dict(
                COMPRESSED_PATH=compressed_path,
                EXTRACT_PATH=extract_path,
                THREADS=config.unpack.sge_slots,
            )
        case _:
            raise ValueError(f"Unknown unpack method: {method}")

    sge.submit(
        str(Path(__file__).parent / "scripts" / f"{method}.sh"),
        env=env,
        queue=config.unpack.sge_queue,
        pe=config.unpack.sge_pe,
        slots=config.unpack.sge_slots,
        name=f"unpack_{compressed_path.name}",
        stderr=config.logdir / f"{compressed_path.name}.{method}.err",
        stdout=config.logdir / f"{compressed_path.name}.{method}.out",
        cwd=compressed_path.parent,
        check=True,
    )


def _extract_callback(
    future: Future,
    /,
    logger: LoggerAdapter,
    samples: data.Samples,
    compressed_path: Path,
    extract_path: Path,
    s_idx: int,
    f_idx: int,
):
    if (exception := future.exception()) is not None:
        logger.error(f"Failed to extract {compressed_path} ({exception})")
        samples[s_idx].files[f_idx] = None
    elif extract_path.exists():
        logger.debug(f"Extracted {extract_path}")
        samples[s_idx].files[f_idx] = extract_path
    elif ((fq1 := extract_path.with_suffix(".gz.1")).exists()) and (
        (fq2 := extract_path.with_suffix(".gz.2")).exists()
    ):
        fq1.rename(fq1.parent / f"{fq1.name.partition('.')[0]}_1.fastq.gz")
        fq2.rename(fq2.parent / f"{fq2.name.partition('.')[0]}_2.fastq.gz")
        logger.debug(f"Extracted {fq1} and {fq2}")
        samples[s_idx].files = [fq1, fq2]
    else:
        logger.error(f"Extraction completed, but {extract_path} does not exist")
        samples[s_idx].files[f_idx] = None


@modules.pre_hook(label="unpack", after=["hcp_fetch"])
def unpack(
    samples: data.Samples,
    config: cfg.Config,
    logger: LoggerAdapter,
    **_,
) -> data.Samples:
    """Extract petagene fasterq files."""
    with ProcessPoolExecutor(config.unpack.parallel) as pool:
        for s_idx, sample in enumerate(samples):
            for f_idx, fastq in enumerate(sample.files):
                if fastq and (compressed_path := Path(fastq)).exists():
                    if compressed_path.suffix == ".fasterq":
                        method = "petagene"
                    elif compressed_path.suffix == ".spring":
                        method = "spring"
                    else:
                        continue

                    extract_path = compressed_path.with_suffix(".fastq.gz")

                    alt_paths = [
                        extract_path.parent
                        / f"{extract_path.name.partition('.')[0]}_1.fastq.gz",
                        extract_path.parent
                        / f"{extract_path.name.partition('.')[0]}_2.fastq.gz",
                    ]

                    if extract_path.exists():
                        logger.debug(f"Extracted file found for {sample.id}")
                        sample.files[f_idx] = extract_path
                        continue
                    elif all(p.exists() for p in alt_paths):
                        logger.debug(f"Extracted files found for {sample.id}")
                        sample.files = alt_paths
                        continue
                    else:
                        logger.info(f"Extracting {compressed_path}")
                        pool.submit(
                            _extract,
                            method,
                            compressed_path=compressed_path,
                            extract_path=extract_path,
                            config=config,
                        ).add_done_callback(
                            partial(
                                _extract_callback,
                                logger=logger,
                                samples=samples,
                                compressed_path=compressed_path,
                                extract_path=extract_path,
                                s_idx=s_idx,
                                f_idx=f_idx,
                            )
                        )

    return samples.__class__([s for s in samples if s is not None])
