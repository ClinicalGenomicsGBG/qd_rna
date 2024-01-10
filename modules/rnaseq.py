"""Module for running nf-core/rnaseq."""

from logging import LoggerAdapter
from pathlib import Path
from cellophane import output, runner, Executor, Config, Samples
from modules.nextflow import nextflow

@output(
    "star_salmon/salmon.*",
    dst_dir="{sample.id}/salmon",
)
@output(
    "star_salmon/{sample.id}",
    dst_dir="{sample.id}/salmon/{sample.id}",
)
@output(
    "star_salmon/{config.rnaseq.aligner}/stringtie",
    dst_dir="{sample.id}/stringtie",
)
@output(
    "multiqc/{config.rnaseq.aligner}",
    dst_dir="{sample.id}/multiqc",
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

    result, _ = nextflow(
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

    result.get()

    return samples
