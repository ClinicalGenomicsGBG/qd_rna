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


def _validate_inputs(config: Config, logger: LoggerAdapter) -> None:
    if any({"genome", x} <= {*config.rnaseq} for x in ["fasta", "gtf", "gene_bed"]):
        logger.warning("Both genome and fasta/gtf/gene_bed provided. Using genome.")

    if (
        "genome" not in config.rnaseq
        and (
            config.rnaseq.aligner == "star_salmon"
            and any(
                not path.exists()
                for path in [
                    config.rnaseq.fasta,
                    config.rnaseq.gtf,
                    config.rnaseq.gene_bed,
                ]
            )
        )
        or (
            config.rnaseq.aligner == "star_rsem"
            and not config.rnaseq.rsem_index.exists()
        )
    ):
        logger.error("Missing required reference files for nf-core/rnaseq.")
        raise SystemExit(1)


def _pipeline_args(
    config: Config,
    nf_samples: Path,
    workdir: Path,
) -> list[str]:
    return [
        f"--outdir {workdir}",
        f"--input {nf_samples}",
        f"--aligner {config.rnaseq.aligner}",
        "--pseudo_aligner salmon",
        (
            f"--salmon_index {config.rnaseq.salmon_index}"
            if "salmon_index" in config.rnaseq
            else ""
        ),
        (
            f"--extra_star_args='--limitSjdbInsertNsj {config.limitSjdbInsertNsj}'"
            if config.rnaseq.aligner == "star_salmon"
            else ""
        ),
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
    ]


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
    logger: LoggerAdapter,
    workdir: Path,
    executor: Executor,
    root: Path,
    **_,
) -> Samples:
    """Run nf-core/rnaseq."""
    if config.rnaseq.skip:
        if not config.copy_skipped:
            samples.output = set()
        return samples

    _validate_inputs(config, logger)
    logger.info("Running nf-core/rnaseq")

    for id_, group in samples.split(by="id"):
        sample_sheet = group.nfcore_samplesheet(
            location=workdir / id_,
            strandedness=config.strandedness,
        )

        nextflow(
            root / "dependencies" / "nf-core" / "rnaseq" / "main.nf",
            *_pipeline_args(
                config=config,
                workdir=workdir / id_,
                nf_samples=sample_sheet,
            ),
            config=config,
            name="rnaseq",
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
            ),
        )

    executor.wait()

    return samples
