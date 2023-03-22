"""Module for running nf-core/rnaseq."""

from logging import LoggerAdapter
from pathlib import Path
from cellophane import cfg, modules, data
from modules.nextflow import nextflow
from modules.qd_rna import Output


@modules.runner(individual_samples=True)
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

    if "rnaseq" in config and not config.rnaseq.skip:
        if any({"genome", x} <= {*config.rnaseq} for x in ["fasta", "gtf", "gene_bed"]):
            logger.warning("Both genome and fasta/gtf/gene_bed provided. Using genome.")

        logger.info("Running nf-core/rnaseq")

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

        src_base = outdir / samples[0].id
        dest_base = Path(samples[0].id)
        samples[0].output = [
            Output(
                src=(src_base / "star_salmon").glob("salmon.*"),
                dest_dir=dest_base / "salmon",
            ),
            Output(
                src=(src_base / "star_salmon" / samples[0].id).glob("*"),
                dest_dir=dest_base / "salmon" / samples[0].id,
            ),
            Output(
                src=(src_base / config.rnaseq.aligner / "stringtie").glob("*"),
                dest_dir=dest_base / "stringtie",
            ),
            Output(
                src=(src_base / "multiqc" / config.rnaseq.aligner).glob("*"),
                dest_dir=dest_base / "multiqc",
            ),
        ]

    return samples