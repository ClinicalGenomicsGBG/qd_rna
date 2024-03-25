"""Module for running nf-core/rnaseq."""

from functools import partial
from logging import LoggerAdapter
from pathlib import Path

from cellophane import Config, Executor, OutputGlob, Samples, output, runner

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
            f"--transcript_fasta {config.rnaseq.transcript_fasta}"
            if "transcript_fasta" in config.rnaseq
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

def _add_optional_outputs(samples: Samples, config: Config) -> None:
    # FIXME: Remove dst_name=None when cellophane is updated from 1.0.0
    if config.rnaseq.aligner == "star_salmon":
        samples.output |= {
            OutputGlob(
                src="{sample.id}/star_salmon/salmon.merged.*",
                dst_dir="{sample.id}/expression/salmon/",
                dst_name=None,
            ),
            OutputGlob(
                src="{sample.id}/star_salmon/{sample.id}",
                dst_dir="{sample.id}/expression/salmon/",
                dst_name=None,
            )
        }

    if config.rnaseq.aligner == "star_rsem":
        samples.output |= {
            OutputGlob(
                src="{sample.id}/star_salmon/rsem.merged.*",
                dst_dir="{sample.id}/expression/rsem/",
                dst_name=None,
            ),
            OutputGlob(
                src="{sample.id}/star_salmon/{sample.id}.genes.results",
                dst_dir="{sample.id}/expression/rsem/",
                dst_name=None,
            ),
            OutputGlob(
                src="{sample.id}/star_salmon/{sample.id}.isoforms.results",
                dst_dir="{sample.id}/expression/rsem/",
                dst_name=None,
            ),
            OutputGlob(
                src="{sample.id}/star_salmon/{sample.id}.stat",
                dst_dir="{sample.id}/expression/rsem/",
                dst_name=None,
            )
        }


@output(
    "{sample.id}/{config.rnaseq.aligner}/{sample.id}.markdup.sorted.bam",
    dst_dir="{sample.id}/expression",
)
@output(
    "{sample.id}/{config.rnaseq.aligner}/{sample.id}.markdup.sorted.bam.bai",
    dst_dir="{sample.id}/expression",
)
@output(
    "{sample.id}/salmon",
    dst_name="{sample.id}/expression/salmon_pseudo",
)
@output(
    "{sample.id}/{config.rnaseq.aligner}/stringtie",
    dst_dir="{sample.id}/expression/",
)
@output(
    "{sample.id}/{config.rnaseq.aligner}/bigwig",
    dst_dir="{sample.id}/expression/",
)
@output(
    "{sample.id}/{config.rnaseq.aligner}/samtools_stats",
    dst_dir="{sample.id}/expression/qc/",
)
@output(
    "{sample.id}/{config.rnaseq.aligner}/picard_metrics",
    dst_dir="{sample.id}/expression/qc/",
)
@output(
    "{sample.id}/{config.rnaseq.aligner}/rseqc",
    dst_dir="{sample.id}/expression/qc/",
)
@output(
    "{sample.id}/{config.rnaseq.aligner}/qualimap",
    dst_dir="{sample.id}/expression/qc/",
)
@output(
    "{sample.id}/{config.rnaseq.aligner}/dupradar",
    dst_dir="{sample.id}/expression/qc/",
)
@output(
    "{sample.id}/{config.rnaseq.aligner}/deseq2_qc",
    dst_dir="{sample.id}/expression/qc/",
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

    _add_optional_outputs(samples, config)

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
            resume=True,
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
