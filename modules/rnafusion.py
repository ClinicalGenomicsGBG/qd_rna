"""Module for running nf-core/rnafusion."""

from logging import LoggerAdapter
from pathlib import Path

from cellophane import Config, Executor, Samples, output, runner

from modules.nextflow import nextflow
from modules.qd_rna import nf_config

# Taken from https://github.com/nf-core/rnafusion/blob/3.0.1/conf/modules.config
# The whole string needs to be included here only to override the limitSjdbInsertNsj
# parameter.
rnafusion_nf_config = """\
process {{
    withName: 'STAR_FOR_STARFUSION' {{
        ext.args = '--twopassMode Basic \
        --outReadsUnmapped None \
        --readFilesCommand zcat \
        --outSAMtype BAM SortedByCoordinate \
        --outSAMstrandField intronMotif \
        --outSAMunmapped Within \
        --chimSegmentMin 12 \
        --chimJunctionOverhangMin 8 \
        --chimOutJunctionFormat 1 \
        --alignSJDBoverhangMin 10 \
        --alignMatesGapMax 100000 \
        --alignIntronMax 100000 \
        --alignSJstitchMismatchNmax 5 -1 5 5 \
        --chimMultimapScoreRange 3 \
        --chimScoreJunctionNonGTAG -4 \
        --chimMultimapNmax 20 \
        --chimNonchimScoreDropMin 10 \
        --peOverlapNbasesMin 12 \
        --peOverlapMMp 0.1 \
        --alignInsertionFlush Right \
        --alignSplicedMateMapLminOverLmate 0 \
        --alignSplicedMateMapLmin 30 \
        --chimOutType Junctions \
        --quantMode GeneCounts \
        --limitSjdbInsertNsj {config.limitSjdbInsertNsj}'
    }}
    withName: 'STAR_FOR_ARRIBA' {{
        ext.args = '--readFilesCommand zcat \
        --outSAMtype BAM Unsorted \
        --outSAMunmapped Within \
        --outBAMcompression 0 \
        --outFilterMultimapNmax 50 \
        --peOverlapNbasesMin 10 \
        --alignSplicedMateMapLminOverLmate 0.5 \
        --alignSJstitchMismatchNmax 5 -1 5 5 \
        --chimSegmentMin 10 \
        --chimOutType WithinBAM HardClip \
        --chimJunctionOverhangMin 10 \
        --chimScoreDropMax 30 \
        --chimScoreJunctionNonGTAG 0 \
        --chimScoreSeparation 1 \
        --chimSegmentReadGapMax 3 \
        --chimMultimapNmax 50 \
        --limitSjdbInsertNsj {config.limitSjdbInsertNsj}'
    }}
}}
"""

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

    if not config.rnafusion.genomes_base.exists():
        logger.error("Missing required reference files for nf-core/rnafusion.")
        logger.error(
            "Instructions for downloading reference files can be found at "
            "https://nf-co.re/rnafusion/3.0.1/docs/usage#download-and-build-references"
        )
        logger.error(
            f"Use {root}/dependencies/nf-core/rnafusion/main.nf when downloading "
            "the reference files"
        )
        raise SystemExit(1)

    nf_config(
        template=rnafusion_nf_config,
        location=workdir / "nextflow.config",
        include=config.nextflow.get("config"),
        config=config,
    )

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
        f"--read_length {config.read_length}",
        f"--tools_cutoff {config.rnafusion.tools_cutoff}",
        f"--fusioncatcher_limitSjdbInsertNsj {config.limitSjdbInsertNsj}",
        f"--fusioninspector_limitSjdbInsertNsj {config.limitSjdbInsertNsj}",
        "--all",
        nxf_config=workdir / "nextflow.config",
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
