"""Module for running nf-core/rnavar."""

from logging import LoggerAdapter
from pathlib import Path

from cellophane import Checkpoints, Config, Executor, Samples, output, runner

from modules.nextflow import nextflow


def _reference_files(config: Config) -> dict[str, Path]:
    return {
        "fasta": Path(config.rnavar.fasta),
        "fasta_fai": Path(config.rnavar.fasta_fai),
        "dict": Path(config.rnavar.dict),
        "gtf": Path(config.rnavar.gtf),
        "star_index": Path(config.rnavar.star_index),
    }


def _validate_inputs(config: Config, logger: LoggerAdapter) -> None:
    missing_options = [
        option
        for option in ("fasta", "fasta_fai", "dict", "gtf", "star_index")
        if not getattr(config.rnavar, option, None)
    ]
    if missing_options:
        logger.error(
            "Missing required nf-core/rnavar configuration options: %s",
            ", ".join(missing_options),
        )
        raise SystemExit(1)

    missing = [
        name for name, path in _reference_files(config).items() if not path.exists()
    ]
    if missing:
        logger.error(
            "Missing required reference files for nf-core/rnavar: %s",
            ", ".join(missing),
        )
        raise SystemExit(1)


def _pipeline_args(config: Config, nf_samples: Path, workdir: Path) -> list[str]:
    references = _reference_files(config)
    return [
        f"--outdir {workdir}",
        f"--input {nf_samples}",
        f"--fasta {references['fasta']}",
        f"--fasta_fai {references['fasta_fai']}",
        f"--dict {references['dict']}",
        f"--gtf {references['gtf']}",
        f"--star_index {references['star_index']}",
        f"--read_length {config.read_length}",
        "--skip_tools baserecalibrator",
    ]


@output(
    "variant_calling/{sample.id}/{sample.id}.haplotypecaller.filtered.vcf.gz",
    dst_dir="{sample.id}_{sample.last_run}_{timestamp}/rnavar",
)
@output(
    "variant_calling/{sample.id}/{sample.id}.haplotypecaller.filtered.vcf.gz.tbi",
    dst_dir="{sample.id}_{sample.last_run}_{timestamp}/rnavar",
)
@output(
    "preprocessing/{sample.id}/{sample.id}.md.bam",
    dst_dir="{sample.id}_{sample.last_run}_{timestamp}/rnavar",
)
@output(
    "preprocessing/{sample.id}/{sample.id}.md.bam.bai",
    dst_dir="{sample.id}_{sample.last_run}_{timestamp}/rnavar",
)
@output(
    "reports",
    dst_name="{sample.id}_{sample.last_run}_{timestamp}/rnavar/reports",
)
@output(
    "pipeline_info",
    dst_name="{sample.id}_{sample.last_run}_{timestamp}/pipeline_info/rnavar",
)
@runner(split_by="id")
def rnavar(
    samples: Samples,
    config: Config,
    logger: LoggerAdapter,
    workdir: Path,
    executor: Executor,
    root: Path,
    checkpoints: Checkpoints,
    **_,
) -> Samples:
    """Run nf-core/rnavar."""
    if config.rnavar.skip:
        samples.output = set()
        return samples

    log_tag = samples[0].id if (n := len(samples.unique_ids)) == 1 else f"{n} samples"
    if checkpoints.main.check():
        logger.info(f"Using previous nf-core/rnavar output ({log_tag})")
        return samples

    _validate_inputs(config, logger)
    logger.info(f"Running nf-core/rnavar ({log_tag})")
    sample_sheet = samples.nfcore_samplesheet(
        location=workdir,
        strandedness=config.strandedness,
        logger=logger,
    )

    nextflow(
        root / "dependencies" / "nf-core" / "rnavar" / "main.nf",
        *_pipeline_args(config, sample_sheet, workdir),
        config=config,
        name="rnavar",
        workdir=workdir,
        resume=True,
        executor=executor,
        conda_spec={"nextflow": ">=26.04"}
    )

    logger.debug(f"nf-core/rnavar finished ({log_tag})")
    checkpoints.main.store()
    return samples
