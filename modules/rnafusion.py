"""Module for running nf-core/rnafusion."""

from logging import LoggerAdapter
from pathlib import Path

from cellophane import Config, Executor, Samples, output, runner

from modules.nextflow import nextflow


def _patch_fusionreport(report_path: Path, sample_id: str):
    """
    Patch fusionreport html to keep fusion data separate from the report itself.

    This is done to make it easier to access the report from the directory listing.
    """
    index_name = f"{sample_id}.fusionreport.html"
    patched_index = (
        (report_path / "index.html")
        .read_text()
        .replace(
            "${fusion.replace('--','_')}.html",
            "fusionreport/${fusion.replace('--','_')}.html",
        )
    )

    Path(index_name).write_text(patched_index, encoding="utf-8")

@output(
    "arriba_visualisation/{sample.id}.pdf",
    dst_dir="{sample.id}",
)
@output(
    "arriba/{sample.id}.*",
    dst_dir="{sample.id}/arriba",
)
@output(
    "fusioncatcher/{sample.id}.*",
    dst_dir="{sample.id}/fusioncatcher",
)
@output(
    "starfusion/{sample.id}.*",
    dst_dir="{sample.id}/starfusion",
)
@output(
    "fusionreport/{sample.id}",
    dst_name="{sample.id}/fusionreport",
)
@output(
    "{sample_id}.fusionreport.html",
    dst_dir="{sample.id}",
)
@output(
    "samtools_sort_for_arriba/{sample.id}_sorted.bam",
    dst_dir="{sample.id}",
)
@output(
    "samtools_index_for_arriba/{sample.id}_sorted.bam.bai",
    dst_dir="{sample.id}",
)
@output(
    "kallisto/{sample.id}.*",
    dst_dir="{sample.id}/kallisto",
)
@output(
    "pipeline_info",
    dst_name="{sample.id}/pipeline_info/rnafusion",
)
@runner()
def rnafusion(
    samples: Samples,
    config: Config,
    workdir: Path,
    executor: Executor,
    logger: LoggerAdapter,
    **_,
) -> Samples:
    """Run nf-core/rnafusion."""
    if config.rnafusion.skip:
        if not config.copy_skipped:
            samples.output = set()
        return samples

    logger.info("Running nf-core/rnafusion")

    sample_sheet = samples.nfcore_samplesheet(
        location=workdir,
        strandedness=config.strandedness,
    )

    nextflow(
        config.rnafusion.nf_main,
        f"--outdir {workdir}",
        f"--input {sample_sheet}",
        f"--genomes_base {config.rnafusion.genomes_base}",
        f"--arriba_ref {config.rnafusion.arriba_ref}",
        f"--arriba_ref_blacklist {config.rnafusion.arriba_blacklist}",
        f"--arriba_ref_protein_domain {config.rnafusion.arriba_protein_domain}",
        f"--read_length {config.read_length}",
        f"--tools_cutoff {config.rnafusion.tools_cutoff}",
        "--fusioncatcher_limitSjdbInsertNsj 4000000",
        "--all",
        config=config,
        name="rnafusion",
        workdir=workdir,
        executor=executor,
    )

    logger.debug(f"nf-core/rnafusion finished for {len(samples)} samples")
    logger.info("Patching fusionreport html")
    for id_, group in samples.split(by="id"):
        logger.debug(f"Patching fusionreport for {id_}")
        try:
            _patch_fusionreport(workdir / f"fusionreport/{id_}", id_)
        except Exception as exception:  # pylint: disable=broad-except
            logger.error(f"Failed to patch fusionreport for {id_}: {exception}")
            for sample in group:
                sample.fail(f"Failed to patch fusionreport: {exception}")

    return samples
