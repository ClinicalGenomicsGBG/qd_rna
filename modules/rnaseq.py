"""Module for running nf-core/rnaseq."""

import time
from logging import LoggerAdapter
from pathlib import Path

from cellophane import cfg, slims, modules, sge


def get_output(aligner: str, outdir: Path):
    """Return a dictionary of output files for the rnaseq module."""
    return {
        aligner: [
            output
            for pattern in [
                "*.bam",
                "*.bai",
                "**/*.bigWig",
                "*.results",
                "*.tsv",
                "*.rds",
            ]
            for output in (outdir / aligner).glob(pattern)
        ],
        "stringtie": [
            output
            for pattern in ["*.gtf", "*.txt", "*.ballgown"]
            for output in (outdir / aligner / "stringtie").glob(pattern)
        ],
        "multiqc": [outdir / "multiqc" / "multiqc_report.html"],
    }


@modules.runner(individual_samples=True)
def rnaseq(
    label: str,
    logger: LoggerAdapter,
    samples: slims.Samples,
    config: cfg.Config,
    scripts_path: Path,
) -> None:
    """Run nf-core/rnaseq."""

    timestamp = time.strftime("%y%m%d-%H%M%S")
    outdir = config.outdir / timestamp / label
    outdir.mkdir(parents=True, exist_ok=True)
    output = get_output(config.rnaseq.aligner, outdir)

    if not config.rnaseq.skip and not config.rnaseq.force and any(output.values()):
        logger.error(f"Found existing nf-core/rnaseq output for {samples[0].id}")
        raise SystemExit(1)

    if not config.rnaseq.skip:
        if any({"genome", x} <= {*config.rnaseq} for x in ["fasta", "gtf", "gene_bed"]):
            logger.warning("Both genome and fasta/gtf/gene_bed provided. Using genome.")

        logger.info("Running nf-core/rnaseq")

        return

        sample_sheet = samples.nfcore_samplesheet(
            location=outdir,
            strandedness=config.rnaseq.strandedness,
        )

        sge.submit(
            str(scripts_path / "nextflow.sh"),
            f"-log {outdir / 'logs' / 'rnaseq.log'}",
            (
                f"-config {config.nextflow.config}"
                if "config" in config.nextflow
                else ""
            ),
            f"run {config.rnaseq.nf_main}",
            "-ansi-log false",
            "-resume" if config.nextflow.resume else "",
            f"-work-dir {config.nextflow.workdir or outdir / 'work'}",
            f"-with-report {outdir / 'logs' / 'rnaseq-execution.html'}",
            f"-profile {config.nextflow.profile}",
            f"--outdir {outdir}",
            f"--input {sample_sheet}",
            f"--aligner {config.rnaseq.aligner}",
            (
                "--pseudo_aligner salmon"
                if config.rnaseq.aligner != "star_salmon"
                else ""
            ),
            (
                "--fasta {config.rnaseq.fasta} "
                "--gtf {config.rnaseq.gtf} "
                "--gene_bed {config.rnaseq.gene_bed}"
                if "genome" not in config.rnaseq
                else f"--genome {config.rnaseq.genome}"
            ),
            (
                f"--star_index {config.rnaseq.index}"
                if config.rnaseq.aligner == "star_salmon"
                else f"--rsem_index {config.rnaseq.index}"
                if config.rnaseq.aligner == "star_rsem"
                else ""
            ),
            env={"_MODULES_INIT": config.modules_init},
            queue=config.nextflow.sge_queue,
            pe=config.nextflow.sge_pe,
            slots=config.nextflow.sge_slots,
            check=True,
            name="rnafusion",
            stderr=config.logdir / f"rnaseq.{timestamp}.err",
            stdout=config.logdir / f"rnaseq.{timestamp}.out",
            cwd=outdir,
        )
