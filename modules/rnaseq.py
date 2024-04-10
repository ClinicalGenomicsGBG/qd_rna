"""Module for running nf-core/rnaseq."""

from logging import LoggerAdapter
from pathlib import Path

from cellophane import (
    Checkpoints,
    Config,
    Executor,
    OutputGlob,
    Samples,
    output,
    runner,
)

from modules.nextflow import nextflow


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
                src="star_salmon/salmon.merged.*",
                dst_dir="{sample.id}_{sample.last_run}/expression/salmon/",
                dst_name=None,
            ),
            OutputGlob(
                src="star_salmon/{sample.id}",
                dst_dir="{sample.id}_{sample.last_run}/expression/salmon/",
                dst_name=None,
            ),
        }

    if config.rnaseq.aligner == "star_rsem":
        samples.output |= {
            OutputGlob(
                src="star_salmon/rsem.merged.*",
                dst_dir="{sample.id}_{sample.last_run}/expression/rsem/",
                dst_name=None,
            ),
            OutputGlob(
                src="star_salmon/{sample.id}.genes.results",
                dst_dir="{sample.id}_{sample.last_run}/expression/rsem/",
                dst_name=None,
            ),
            OutputGlob(
                src="star_salmon/{sample.id}.isoforms.results",
                dst_dir="{sample.id}_{sample.last_run}/expression/rsem/",
                dst_name=None,
            ),
            OutputGlob(
                src="star_salmon/{sample.id}.stat",
                dst_dir="{sample.id}_{sample.last_run}/expression/rsem/",
                dst_name=None,
            ),
        }


@output(
    "{config.rnaseq.aligner}/{sample.id}.markdup.sorted.bam",
    dst_dir="{sample.id}_{sample.last_run}/expression",
)
@output(
    "{config.rnaseq.aligner}/{sample.id}.markdup.sorted.bam.bai",
    dst_dir="{sample.id}_{sample.last_run}/expression",
)
@output(
    "salmon",
    dst_name="{sample.id}_{sample.last_run}/expression/salmon_pseudo",
)
@output(
    "{config.rnaseq.aligner}/stringtie",
    dst_dir="{sample.id}_{sample.last_run}/expression/",
)
@output(
    "{config.rnaseq.aligner}/bigwig",
    dst_dir="{sample.id}_{sample.last_run}/expression/",
)
@output(
    "{config.rnaseq.aligner}/samtools_stats",
    dst_dir="{sample.id}_{sample.last_run}/expression/qc/",
)
@output(
    "{config.rnaseq.aligner}/picard_metrics",
    dst_dir="{sample.id}_{sample.last_run}/expression/qc/",
)
@output(
    "{config.rnaseq.aligner}/rseqc",
    dst_dir="{sample.id}_{sample.last_run}/expression/qc/",
)
@output(
    "{config.rnaseq.aligner}/qualimap",
    dst_dir="{sample.id}_{sample.last_run}/expression/qc/",
)
@output(
    "{config.rnaseq.aligner}/dupradar",
    dst_dir="{sample.id}_{sample.last_run}/expression/qc/",
)
@output(
    "{config.rnaseq.aligner}/deseq2_qc",
    dst_dir="{sample.id}_{sample.last_run}/expression/qc/",
)
@output(
    "multiqc/{config.rnaseq.aligner}",
    dst_name="{sample.id}_{sample.last_run}/multiqc",
)
@output(
    "pipeline_info",
    dst_name="{sample.id}_{sample.last_run}/pipeline_info/rnaseq",
)
@runner(split_by="id")
def rnaseq(
    samples: Samples,
    config: Config,
    logger: LoggerAdapter,
    workdir: Path,
    executor: Executor,
    root: Path,
    checkpoints: Checkpoints,
    **_,
) -> Samples:
    """Run nf-core/rnaseq."""
    if config.rnaseq.skip:
        samples.output = set()
        return samples

    _add_optional_outputs(samples, config)

    if checkpoints.main.check():
        logger.info(f"Using previous nf-core/rnaseq output ({samples[0].id})")
        return samples

    _validate_inputs(config, logger)
    logger.info("Running nf-core/rnaseq ({samples[0].id})")

    sample_sheet = samples.nfcore_samplesheet(
        location=workdir,
        strandedness=config.strandedness,
    )

    nextflow(
        root / "dependencies" / "nf-core" / "rnaseq" / "main.nf",
        *_pipeline_args(
            config=config,
            workdir=workdir,
            nf_samples=sample_sheet,
        ),
        config=config,
        name="rnaseq",
        workdir=workdir,
        resume=True,
        executor=executor,
    )

    logger.debug(f"nf-core/rnaseq finished for sample {samples[0].id}")
    checkpoints.main.store()

    return samples
