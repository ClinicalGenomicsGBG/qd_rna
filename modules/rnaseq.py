"""Module for running nf-core/rnaseq."""

from logging import LoggerAdapter
from pathlib import Path
from cellophane import output, runner, Executor, Config, Samples
from modules.nextflow import nextflow

@output(
    "{config.rnaseq.aligner}",
    dst_name="{sample.id}/{config.rnaseq.aligner}",
)
@output(
    "multiqc/{config.rnaseq.aligner}",
    dst_name="{sample.id}/multiqc",
)
@output(
    "pipeline_info",
    dst_name="{sample.id}/pipeline_info/rnaseq",
)
@runner(individual_samples=True, link_by="id")
def rnaseq(
    samples: Samples,
    config: Config,
    label: str,
    logger: LoggerAdapter,
    workdir: Path,
    executor: Executor,
    **_,
) -> Samples:
    """Run nf-core/rnaseq."""

    if config.rnaseq.skip:
        if not config.copy_skipped:
            samples.output = set()
        return samples
    
    if any({"genome", x} <= {*config.rnaseq} for x in ["fasta", "gtf", "gene_bed"]):
        logger.warning("Both genome and fasta/gtf/gene_bed provided. Using genome.")

    sample_sheet = samples.nfcore_samplesheet(
        location=workdir,
        strandedness=config.strandedness,
    )

    nextflow(
        config.rnaseq.nf_main,
        f"--outdir {workdir}",
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
        name=label,
        workdir=workdir,
        executor=executor,
    )

    return samples
