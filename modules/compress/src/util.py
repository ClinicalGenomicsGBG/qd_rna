from __future__ import annotations
from fnmatch import fnmatch

from pathlib import Path
from logging import LoggerAdapter
import shlex

from cellophane import Executor


def should_convert(out, rules, workdir) -> tuple[bool, str | None, Path | None]:
    src = Path(out.src)
    src_key = str(src.relative_to(workdir))

    for rule in rules:
        if fnmatch(src_key, str(rule.pattern)):
            return True, rule.method.lower(), Path(rule.reference)
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

    if name is None:
        name = f"cram_{input_bam.stem}"

    if conda_spec is None:
        conda_spec = {
            "channels": ["bioconda", "conda-forge"],
            "dependencies": ["samtools=1.16"],
        }

    cmd = (
        "samtools view "
        f"-@ {int(threads)} "
        f"-C -T {shlex.quote(str(reference))} "
        f"-o {shlex.quote(str(output_cram))} "
        f"{shlex.quote(str(input_bam))}"
        f" && samtools index -@ {int(threads)} {shlex.quote(str(output_cram))}"
    )

    if logger:
        logger.info(f"Submitting CRAM compression: {input_bam.name} -> {output_cram.name}")

    return executor.submit(
        "bash",
        "-lc",
        cmd,
        name=name,
        workdir=workdir,
        cpus=int(threads),
        conda_spec=conda_spec,
        callback=callback,
        error_callback=error_callback,
    )

