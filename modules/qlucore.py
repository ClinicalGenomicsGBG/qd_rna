"""Module for running nf-core/rnaseq."""

from logging import LoggerAdapter
from pathlib import Path
from copy import deepcopy

from cellophane import cfg, modules
from modules.slims import SlimsSamples

from modules.nextflow import nextflow


@modules.runner()
def qlucore(
    samples: SlimsSamples,
    config: cfg.Config,
    timestamp: str,
    label: str,
    logger: LoggerAdapter,
    root: Path,
    outdir: Path,
) -> None:
    """Run nf-core/rnaseq (Mapping for qlucore)."""

    _samples = deepcopy(samples)
    for idx, s in enumerate(samples):
        _samples[idx].id = f"{s.id}_{s.run}" if "run" in s and s.run else s.id

    if "qlucore" in config and not config.qlucore.skip:
        if any(
            {"genome", x} <= {*config.qlucore} for x in ["fasta", "gtf", "gene_bed"]
        ):
            logger.warning("Both genome and fasta/gtf/gene_bed provided. Using genome.")

        logger.info("Running nf-core/rnaseq (for qlucore)")

        sample_sheet = samples.nfcore_samplesheet(
            location=outdir,
            strandedness=config.rnaseq.strandedness,
        )

        if "workdir" in config.nextflow:
            config.nextflow.workdir.mkdir(parents=True, exist_ok=True)

        nextflow(
            config.rnaseq.nf_main,
            f"--outdir {outdir}",
            f"--input {sample_sheet}",
            "--aligner star_salmon",
            "--skip_qc",
            "--skip_bigwig",
            f"--star_index {config.qlucore.star_index}",
            (
                f"--fasta {config.qlucore.fasta} "
                f"--gtf {config.qlucore.gtf} "
                f"--gene_bed {config.qlucore.gene_bed}"
                if "genome" not in config.qlucore
                else f"--genome {config.qlucore.genome}"
            ),
            config=config,
            name="qlucore",
            workdir=outdir / "work",
            stderr=outdir / "logs" / f"{label}.{timestamp}.err",
            stdout=outdir / "logs" / f"{label}.{timestamp}.out",
            cwd=outdir,
        )
