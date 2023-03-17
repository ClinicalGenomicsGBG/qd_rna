"""Module for running nf-core/rnaseq."""

from logging import LoggerAdapter
from pathlib import Path
from copy import deepcopy

from cellophane import cfg, modules
from modules.slims import SlimsSamples
from modules.nextflow import nextflow
from modules.qd_rna import Output


def get_output(aligner: str, outdir: Path):
    """Return a dictionary of output files for the rnaseq module."""
    return {
        aligner: [
            output
            for pattern in [
                "*.bam",
                "*.bai",
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
    samples: SlimsSamples,
    config: cfg.Config,
    timestamp: str,
    label: str,
    logger: LoggerAdapter,
    root: Path,
    outdir: Path,
) -> None:
    """Run nf-core/rnaseq."""

    _sample = deepcopy(samples)
    _sample[0].id = (
        f"{_sample[0].id}_{_sample[0].run}"
        if "run" in _sample[0] and _sample[0].run
        else _sample[0].id
    )

    if "rnaseq" in config and not config.rnaseq.skip:
        if any({"genome", x} <= {*config.rnaseq} for x in ["fasta", "gtf", "gene_bed"]):
            logger.warning("Both genome and fasta/gtf/gene_bed provided. Using genome.")

        logger.info("Running nf-core/rnaseq")

        sample_sheet = _sample.nfcore_samplesheet(
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
        samples.output += [
            Output(
                src = (outdir / _sample[0].id / config.rnaseq.aligner).glob(f"{_sample[0].id}.markdup.sorted.bam*"),
                dest = Path(_sample[0].id),
            ),
            Output(
                src = (outdir / _sample[0].id / "star_salmon").glob("salmon.*"),
                dest = Path(_sample[0].id) / "salmon",
            ),
            Output(
                src = (outdir / _sample[0].id / "star_salmon" / _sample[0].id).glob("*"),
                dest = Path(_sample[0].id) / "salmon" / _sample[0].id,
            ),
            Output(
                src = (outdir / _sample[0].id / config.rnaseq.aligner / "stringtie").glob("*"),
                dest = Path(_sample[0].id) / "stringtie",
            ),
            Output(
                src = (outdir / _sample[0].id / "multiqc" / config.rnaseq.aligner).glob("*"),
                dest = Path(_sample[0].id) / "multiqc",
            )
        ]

        return samples
    