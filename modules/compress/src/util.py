from __future__ import annotations
from fnmatch import fnmatch

from pathlib import Path
from logging import LoggerAdapter

from cellophane import Executor

ROOT = Path(__file__).parent.parent


def should_convert(out, rules, workdir, logger) -> tuple[bool, str | None, Path | None]:
    src = Path(out.src)
    src_key = str(src.relative_to(workdir))

    for rule in rules:
        pattern = rule.get("pattern")
        if not pattern:
            logger.warning(f"Compression rule missing pattern: {rule}")
            continue

        if fnmatch(src_key, str(pattern)):
            logger.debug(f"Compression rule matched for '{src_key}': {rule}")
            method = rule.get("method")
            if not method:
                raise ValueError(f"Compression rule matched '{src_key}' but no method was provided")

            method = str(method).lower()

            if method == "cram":
                ref = rule.get("reference")
                if not ref:
                    raise ValueError(f"CRAM rule matched '{src_key}' but no reference was provided")
                return True, method, Path(str(ref))

            raise ValueError(f"Unsupported compression method '{method}' for matched rule on '{src_key}'")

    return False, None, None


def find_bai_output(outputs: set, bam_src: Path):
    bai1 = bam_src.with_suffix(".bam.bai")
    bai2 = bam_src.with_suffix(".bai")
    for o in outputs:
        s = Path(o.src)
        if s == bai1 or s == bai2:
            return o
    return None


def cram_compress(
    *,
    executor: Executor,
    input_bam: str | Path,
    reference: str | Path,
    threads: int = 4,
    logger: LoggerAdapter | None = None,
    output_cram: str | Path | None = None,
    name: str | None = None,
    workdir: Path | None = None,
    conda_spec: dict | None = None,
    callback=None,
    error_callback=None,
):
    input_bam = Path(input_bam)
    reference = Path(reference)
    output_cram = Path(output_cram) if output_cram else input_bam.with_suffix(".cram")
    threads_i = int(threads)

    if name is None:
        name = f"cram_{input_bam.stem}"

    if conda_spec is None:
        conda_spec = {
            "channels": ["bioconda", "conda-forge"],
            "dependencies": ["samtools=1.16"],
        }

    script = str(ROOT / "scripts" / "cram_compress.sh")

    if logger:
        logger.info(f"Submitting CRAM compression: {input_bam.name} -> {output_cram.name}")
        logger.debug(
            "CRAM script command: "
            f"{script} {input_bam} {output_cram} {reference} {threads_i}"
        )

    return executor.submit(
        str(script),
        str(input_bam),
        str(output_cram),
        str(reference),
        str(threads_i),
        name=name,
        workdir=workdir,
        cpus=threads_i,
        # conda_spec=conda_spec,  # Use module samtools in the script instead
        callback=callback,
        error_callback=error_callback,
    )
