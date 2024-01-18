"""Module for running nf-core/rnaseq."""

from functools import partial
from logging import LoggerAdapter
from pathlib import Path

from cellophane import Config, Executor, Samples, output, runner

from modules.nextflow import nextflow


def _callback(
    result: None,
    /,
    sample_id: str,
    logger: LoggerAdapter,
):
    del result  # unused
    logger.debug(f"nf-core/rnaseq finished for sample {sample_id}")


def _error_callback(
    exception: Exception,
    /,
    logger: LoggerAdapter,
    sample_id: str,
    group: Samples,
):
    reason = f"nf-core/rnaseq failed for {sample_id}: {exception}"
    logger.error(reason)
    for sample in group:
        sample.fail(reason)

@output(
    "{sample.id}/{config.rnaseq.aligner}",
    dst_name="{sample.id}/{config.rnaseq.aligner}",
)
@output(
    "{sample.id}/multiqc/{config.rnaseq.aligner}",
    dst_name="{sample.id}/multiqc",
)
@output(
    "{sample.id}/pipeline_info",
    dst_name="{sample.id}/pipeline_info/rnaseq",
)
@runner()
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

    logger.info("Running nf-core/rnaseq")

    if any({"genome", x} <= {*config.rnaseq} for x in ["fasta", "gtf", "gene_bed"]):
        logger.warning("Both genome and fasta/gtf/gene_bed provided. Using genome.")

    for id_, group in samples.split(link_by="id"):
        sample_sheet = group.nfcore_samplesheet(
            location=workdir / id_,
            strandedness=config.strandedness,
        )

        nextflow(
            config.rnaseq.nf_main,
            f"--outdir {workdir / id_}",
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
            workdir=workdir / id_,
            executor=executor,
            check=False,
            callback=partial(
                _callback,
                sample_id=id_,
                logger=logger,
            ),
            error_callback=partial(
                _error_callback,
                group=group,
                sample_id=id_,
                logger=logger,
            )
        )

    executor.wait()

    return samples
