"""Module for running nf-core/rnaseq."""

from logging import LoggerAdapter
from pathlib import Path
from cellophane import cfg, modules, data
from modules.nextflow import nextflow


@modules.runner(individual_samples=True, link_by="id")
def rnaseq(
    samples: data.Samples,
    config: cfg.Config,
    timestamp: str,
    label: str,
    logger: LoggerAdapter,
    root: Path,
    outdir: Path,
) -> data.Samples:
    """Run nf-core/rnaseq."""

    if not config.rnaseq.skip:
        if any({"genome", x} <= {*config.rnaseq} for x in ["fasta", "gtf", "gene_bed"]):
            logger.warning("Both genome and fasta/gtf/gene_bed provided. Using genome.")

        sample_sheet = samples.nfcore_samplesheet(
            location=outdir,
            strandedness=config.strandedness,
        )

        if "workdir" in config.nextflow:
            config.nextflow.workdir.mkdir(parents=True, exist_ok=True)

        nextflow(
            config.rnaseq.nf_main,
            f"--outdir {outdir}",
            f"--input {sample_sheet}",
            f"--aligner {config.rnaseq.aligner}",
            "--pseudo_aligner salmon",
            (
                f"--fasta {config.rnaseq.fasta} "
                f"--gtf {config.rnaseq.gtf} "
                f"--gene_bed {config.rnaseq.gene_bed}"
                if "genome" not in config.rnaseq
                else f"--genome {config.rnaseq.genome}"
            ),
            (
                f"--rsem_index {config.rnaseq.rsem_index}"
                if config.rnaseq.aligner == "star_rsem"
                else f"--star_index {config.rnaseq.star_index}"
            ),
            config=config,
            name="rnaseq",
            workdir=outdir / "work",
            report=outdir / "logs" / f"{label}.{timestamp}.nextflow_report.html",
            log=outdir / "logs" / f"{label}.{timestamp}.nextflow.log",
            stderr=outdir / "logs" / f"{label}.{timestamp}.nextflow.err",
            stdout=outdir / "logs" / f"{label}.{timestamp}.nextflow.out",
            cwd=outdir,
        )

    if not config.rnaseq.skip or config.results.copy_skipped:
        for sample in samples:
            sample.output = [
                data.Output(
                    src=(outdir / "star_salmon").glob("salmon.*"),
                    dest_dir=Path(samples[0].id) / "salmon",
                ),
                data.Output(
                    src=outdir / "star_salmon" / samples[0].id,
                    dest_dir=Path(samples[0].id) / "salmon" / samples[0].id,
                ),
                data.Output(
                    src=outdir / config.rnaseq.aligner / "stringtie",
                    dest_dir=Path(samples[0].id) / "stringtie",
                ),
                data.Output(
                    src=outdir / "multiqc" / config.rnaseq.aligner,
                    dest_dir=Path(samples[0].id) / "multiqc",
                ),
            ]

    return samples
