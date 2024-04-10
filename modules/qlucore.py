"""Module for running nf-core/rnaseq."""

import re
from functools import partial
from logging import LoggerAdapter
from pathlib import Path

from cellophane import (
    Checkpoint,
    Checkpoints,
    Config,
    Executor,
    Samples,
    output,
    runner,
)

from modules.nextflow import nextflow
from modules.qd_rna import nf_config

qlucore_data = """\
<Header Producer='Qlucore' Format='PatientData' FormatVersion='0.1' QFFVersion='1.0'/>
<PatientData>
  <PatientName>PATIENT NAME</PatientName>
  <PatientId>{id}</PatientId>
  <SampleOrigin>{id}</SampleOrigin>
  <SampleTissue>Blood sample</SampleTissue>
  <Technology>RNA Seq.</Technology>
</PatientData>
"""

qlucore_nf_config = """
process {{
  withName: 'STAR_FOR_STARFUSION' {{
    ext.args = [
        '--twopassMode Basic',
        '--outReadsUnmapped None',
        '--readFilesCommand zcat',
        '--outSAMtype BAM SortedByCoordinate',
        '--outSAMstrandField intronMotif',
        '--outSAMunmapped Within',
        '--chimSegmentMin 12',
        '--chimJunctionOverhangMin 8',
        '--chimOutJunctionFormat 1',
        '--alignSJDBoverhangMin 10',
        '--alignMatesGapMax 100000',
        '--alignIntronMax 100000',
        '--alignSJstitchMismatchNmax 5 -1 5 5',
        '--chimMultimapScoreRange 3',
        '--chimScoreJunctionNonGTAG -4',
        '--chimMultimapNmax 20',
        '--chimNonchimScoreDropMin 10',
        '--peOverlapNbasesMin 12',
        '--peOverlapMMp 0.1',
        '--alignInsertionFlush Right',
        '--alignSplicedMateMapLminOverLmate 0',
        '--alignSplicedMateMapLmin 30',
        '--chimOutType Junctions',
        '--outFilterMultimapNmax 200',
        '--limitSjdbInsertNsj {config.limitSjdbInsertNsj}'
    ].join(' ').trim()
  }}
}}
"""


def _subsample_callback(
    result: None,
    /,
    logger: LoggerAdapter,
    sample_id: str,
    checkpoint: Checkpoint,
):
    del result  # unused
    logger.debug(f"Subsampling finished for {sample_id}")
    checkpoint.store(qlucore_nf_config)


def _subsample_error_callback(
    exception: Exception,
    /,
    logger: LoggerAdapter,
    sample_id: str,
    group: Samples,
):
    reason = f"Subsampling failed for {sample_id} - {exception}"
    logger.error(reason)
    for sample in group:
        sample.fail(reason)


def _calculate_subsample_frac(
    id_: str,
    workdir: Path,
    config: Config,
    logger: LoggerAdapter,
) -> float:
    try:
        star_log = workdir / "star_for_starfusion" / f"{id_}.Log.final.out"

        if match_ := re.search(
            r"Uniquely mapped reads number[^0-9]+([0-9]+)",
            star_log.read_text(),
        ):
            mapped = int(match_.groups()[0])
    except Exception as exc:  # pylint: disable=broad-except
        logger.warning(
            "Unhandled exception in subsample fraction calculation",
            exc_info=exc,
        )
        return config.qlucore.subsample.fallback_fraction

    return config.qlucore.subsample.target / mapped


def _subsample(
    samples: Samples,
    config: Config,
    logger: LoggerAdapter,
    root: Path,
    workdir: Path,
    executor: Executor,
    checkpoints: Checkpoints,
):
    logger.info("Subsampling output BAM(s)")
    for id_, group in samples.split(by="id"):
        output_path = Path(
            workdir
            / "star_for_starfusion"
            / f"{id_}.Aligned.sortedByCoord.out.downsampled.bam"
        )

        if output_path.exists() and checkpoints[f"downsample_{id_}"].check(
            qlucore_nf_config
        ):
            logger.info(f"Using previous downsampled BAM for {id_}")
            continue

        frac = _calculate_subsample_frac(
            id_=id_,
            workdir=workdir,
            config=config,
            logger=logger,
        )

        if frac >= 1.0:
            logger.warning(f"Sample {id_} has too few reads to subsample")
            continue

        logger.debug(f"Subsampling {id_} with fraction {frac}")
        executor.submit(
            str(root / "scripts" / "qlucore_subsample.sh"),
            name=f"qlucore_subsample_{id_}",
            workdir=workdir / id_,
            cpus=config.qlucore.subsample.threads,
            env={
                "_QLUCORE_SUBSAMPLE_INIT": config.qlucore.subsample.init,
                "_QLUCORE_SUBSAMPLE_THREADS": config.qlucore.subsample.threads,
                "_QLUCORE_SUBSAMPLE_FRAC": frac,
                "_QLUCORE_SUBSAMPLE_INPUT_BAM": (
                    workdir
                    / "star_for_starfusion"
                    / f"{id_}.Aligned.sortedByCoord.out.bam"
                ),
            },
            callback=partial(
                _subsample_callback,
                logger=logger,
                sample_id=id_,
                checkpoint=checkpoints[f"downsample_{id_}"],
            ),
            error_callback=partial(
                _subsample_error_callback,
                logger=logger,
                sample_id=id_,
                group=group,
            ),
        )

    executor.wait()


def _validate_inputs(
    config: Config,
    workdir: Path,
    logger: LoggerAdapter,
) -> None:
    (workdir / "dummy.fa").touch()
    (workdir / "dummy.fai").touch()
    (workdir / "dummy.gtf").touch()
    (workdir / "dummy.refflat").touch()
    (workdir / "dummy.interval_list").touch()

    if any(
        not path.exists()
        for path in [
            config.qlucore.starfusion_ref,
            config.qlucore.starfusion_ref / "ref_genome.fa.star.idx",
        ]
    ):
        logger.error("STARFusion reference not found")
        raise SystemExit(1)


def _pipeline_args(
    config: Config,
    workdir: Path,
    nf_samples: Path,
):
    return [
        "--starfusion",
        "--skip_qc",
        "--skip_vis",
        "--star_ignore_sjdbgtf",
        f"--fasta {workdir / 'dummy.fa'}",
        f"--fai {workdir / 'dummy.fai'}",
        f"--transcript {workdir / 'dummy.fa'}",
        f"--gtf {workdir / 'dummy.gtf'}",
        f"--chrgtf {workdir / 'dummy.gtf'}",
        f"--refflat {workdir / 'dummy.refflat'}",
        f"--rrna_intervals {workdir / 'dummy.interval_list'}",
        f"--outdir {workdir}",
        f"--input {nf_samples}",
        f"--starfusion_ref {config.qlucore.starfusion_ref}",
        f"--starindex_ref {config.qlucore.starfusion_ref}/ref_genome.fa.star.idx",
        f"--read_length {config.read_length}",
    ]


@output(
    "starfusion/{sample.id}.*.tsv",
    dst_dir="{sample.id}_{sample.last_run}/qlucore",
)
@output(
    "star_for_starfusion/{sample.id}.Aligned.sortedByCoord.out.bam",
    dst_dir="{sample.id}_{sample.last_run}/qlucore",
)
@output(
    "star_for_starfusion/{sample.id}.Aligned.sortedByCoord.out.bam.bai",
    dst_dir="{sample.id}_{sample.last_run}/qlucore",
)
@output(
    "star_for_starfusion/{sample.id}.Aligned.sortedByCoord.out.downsampled.bam",
    dst_dir="{sample.id}_{sample.last_run}/qlucore",
    checkpoint="downsample_{sample.id}",
    optional=True,
)
@output(
    "{sample.id}.qlucore.txt",
    dst_dir="{sample.id}_{sample.last_run}/qlucore",
    checkpoint="qlucore",
)
@runner()
def qlucore(
    samples: Samples,
    config: Config,
    logger: LoggerAdapter,
    root: Path,
    workdir: Path,
    executor: Executor,
    checkpoints: Checkpoints,
    **_,
) -> None:
    """Run nf-core/rnaseq (Mapping for qlucore)."""

    if config.qlucore.skip:
        if not config.copy_skipped:
            samples.output = set()
        return samples

    if checkpoints.main.check(qlucore_nf_config):
        logger.info("Using previous nf-core/rnafusion output")
    else:
        _validate_inputs(
            config=config,
            workdir=workdir,
            logger=logger,
        )

        sample_sheet = samples.nfcore_samplesheet(
            location=workdir,
            strandedness=config.strandedness,
        )

        nf_config(
            template=qlucore_nf_config,
            location=workdir / "nextflow.config",
            include=config.nextflow.get("config"),
            config=config,
        )

        logger.info("Running nf-core/rnafusion")

        nextflow(
            root / "dependencies" / "nf-core" / "rnafusion" / "main.nf",
            *_pipeline_args(
                config=config,
                workdir=workdir,
                nf_samples=sample_sheet,
            ),
            nxf_config=workdir / "nextflow.config",
            config=config,
            name="qlucore",
            workdir=workdir,
            resume=True,
            executor=executor,
        )

        logger.debug(f"nf-core/rnafusion finished for {len(samples)} samples")
        checkpoints.main.store(qlucore_nf_config)

    _subsample(
        samples=samples,
        config=config,
        logger=logger,
        root=root,
        workdir=workdir,
        executor=executor,
        checkpoints=checkpoints,
    )

    if checkpoints.qlucore.check(qlucore_data):
        logger.info("Using previous qlucore output")
    else:
        for id_, _ in samples.split(by="id"):
            (workdir / f"{id_}.qlucore.txt").write_text(qlucore_data.format(id=id_))
        checkpoints.qlucore.store(qlucore_data)

    return samples
