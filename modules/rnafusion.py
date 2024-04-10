"""Module for running nf-core/rnafusion."""

from logging import LoggerAdapter
from pathlib import Path

from cellophane import Checkpoints, Config, Executor, Samples, output, runner

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


def _patch_fusionreport(samples: Samples, workdir: Path, logger: LoggerAdapter):
    """
    Patch fusionreport html to keep fusion data separate from the report itself.

    This is done to make it easier to access the report from the directory listing.
    """
    logger.info("Patching fusionreport html")
    for id_, group in samples.split(by="id"):
        logger.debug(f"Patching fusionreport for {id_}")
        original = workdir / f"fusionreport/{id_}/{id_}_fusionreport_index.html"
        patched = workdir / f"{id_}.fusionreport.html"

        try:
            index = original.read_text(encoding="utf-8").replace(
                "${fusion.replace('--','_')}.html",
                "fusionreport/${fusion.replace('--','_')}.html",
            )

            patched.write_text(index, encoding="utf-8")
        except Exception as exception:  # pylint: disable=broad-except
            logger.error(f"Failed to patch fusionreport for {id_}: {exception}")
            for sample in group:
                sample.fail(f"Failed to patch fusionreport: {exception}")


def _validate_inputs(config: Config, root: Path, logger: LoggerAdapter):
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


def _pipeline_args(config: Config, workdir: Path, nf_samples: Path, /):
    return [
        f"--outdir {workdir}",
        f"--input {nf_samples}",
        f"--genomes_base {config.rnafusion.genomes_base}",
        f"--read_length {config.read_length}",
        f"--tools_cutoff {config.rnafusion.tools_cutoff}",
        f"--fusioncatcher_limitSjdbInsertNsj {config.limitSjdbInsertNsj}",
        f"--fusioninspector_limitSjdbInsertNsj {config.limitSjdbInsertNsj}",
        "--all",
    ]


@output(
    "arriba_visualisation/{sample.id}_combined_fusions_arriba_visualisation.pdf",
    dst_dir="{sample.id}_{sample.last_run}",
    checkpoint="DUMMY",
)
@output(
    "arriba/{sample.id}.*",
    dst_dir="{sample.id}_{sample.last_run}/arriba",
)
@output(
    "fusioncatcher/{sample.id}.*",
    dst_dir="{sample.id}_{sample.last_run}/fusioncatcher",
)
@output(
    "starfusion/{sample.id}.*",
    dst_dir="{sample.id}_{sample.last_run}/starfusion",
)
@output(
    "fusionreport/{sample.id}",
    dst_name="{sample.id}_{sample.last_run}/fusionreport",
)
@output(
    "{sample.id}.fusionreport.html",
    dst_dir="{sample.id}_{sample.last_run}",
)
@output(
    "star_for_starfusion/{sample.id}.Aligned.sortedByCoord.out.bam",
    dst_dir="{sample.id}_{sample.last_run}",
)
@output(
    "star_for_starfusion/{sample.id}.Aligned.sortedByCoord.out.bam.bai",
    dst_dir="{sample.id}_{sample.last_run}",
)
@output(
    "pipeline_info",
    dst_name="{sample.id}_{sample.last_run}/pipeline_info/rnafusion",
)
@runner()
def rnafusion(
    samples: Samples,
    config: Config,
    workdir: Path,
    executor: Executor,
    logger: LoggerAdapter,
    root: Path,
    checkpoints: Checkpoints,
    **_,
) -> Samples:
    """Run nf-core/rnafusion."""

    if config.rnafusion.skip:
        if not config.copy_skipped:
            samples.output = set()
        return samples

    if checkpoints.main.check(rnafusion_nf_config):
        logger.info("Using previous nf-core/rnafusion output")
        return samples


    _validate_inputs(
        config=config,
        root=root,
        logger=logger,
    )

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
        root / "dependencies" / "nf-core" / "rnafusion" / "main.nf",
        *_pipeline_args(config, workdir, sample_sheet),
        nxf_config=workdir / "nextflow.config",
        config=config,
        name="rnafusion",
        workdir=workdir,
        resume=True,
        executor=executor,
    )

    logger.debug(f"nf-core/rnafusion finished for {len(samples)} samples")

    _patch_fusionreport(
        samples=samples,
        workdir=workdir,
        logger=logger,
    )

    checkpoints.main.store(rnafusion_nf_config)

    return samples
