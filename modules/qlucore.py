"""Module for running nf-core/rnaseq."""

from logging import LoggerAdapter
from pathlib import Path

from cellophane import cfg, modules, data

from modules.qd_rna import Output
from modules.nextflow import nextflow

@modules.runner()
def qlucore(
    samples: data.Samples,
    config: cfg.Config,
    timestamp: str,
    label: str,
    logger: LoggerAdapter,
    root: Path,
    outdir: Path,
) -> data.Samples:
    """Run nf-core/rnaseq (Mapping for qlucore)."""

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
            report=outdir / "logs" / f"{label}.{timestamp}.nextflow_report.html",
            log=outdir / "logs" / f"{label}.{timestamp}.nextflow.log",
            stderr=outdir / "logs" / f"{label}.{timestamp}.nextflow.err",
            stdout=outdir / "logs" / f"{label}.{timestamp}.nextflow.out",
            cwd=outdir,
        )

        for idx, sample in enumerate(samples):
            samples[idx].output = [
                Output(
                    src = (outdir / "star_salmon").glob(f"{sample.id}.*.bam*"),
                    dest_dir = Path(sample.id) / "qlucore",
                ),
            ]
        
    return samples
