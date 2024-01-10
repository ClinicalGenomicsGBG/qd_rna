"""Module for running nf-core/rnaseq."""

from mpire.async_result import AsyncResult
from logging import LoggerAdapter
from pathlib import Path
from functools import partial

from cellophane import output, runner, Executor, Config, Sample, Samples

from modules.nextflow import nextflow


qlucore_data = """\
<Header Producer='Qlucore' Format='PatientData' FormatVersion='0.1' QFFVersion='1.0'/>
<PatientData>
  <PatientName>PATIENT NAME</PatientName>
  <PatientId>{id}</PatientId>
  <SampleOrigin>{meta.run}</SampleOrigin>
  <SampleTissue>Blood sample</SampleTissue>
  <Technology>RNA Seq.</Technology>
</PatientData>
"""

qlucore_nf_config = """\
process {
  withName: 'STAR_FOR_STARFUSION' {
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
        '--limitSjdbInsertNsj 4000000'
    ].join(' ').trim()
  }
}
"""


def _subsample_callback(
    result: AsyncResult,
    /,
    logger: LoggerAdapter,
    workdir: Path,
    sample: Sample,
):
    del result  # unused
    logger.info(f"Subsampling finished for {sample.id}")
    with open(workdir / f"{sample.id}.qlucore.txt", "w") as f:
        f.write(qlucore_data.format(id=sample.id, meta=sample.meta or ""))


def _subsample_error_callback(
    result: AsyncResult,
    /,
    logger: LoggerAdapter,
    sample: Sample,
):
    try:
        result.get()
    except Exception as exception:
        reason = f"Subsampling failed for {sample.id} - {exception}"
    else:
        reason = f"Subsampling failed for {sample.id}"
    logger.error(reason)
    sample.fail(reason)


@output(
    "{sample.id}.qlucore.txt",
    dst_dir="{sample.id}/qlucore",
)
@output(
    "star_for_starfusion/{sample.id}.*.ba*",
    dst_dir="{sample.id}/qlucore",
)
@output(
    "starfusion/{sample.id}.*.tsv",
    dst_dir="{sample.id}/qlucore",
)
@runner()
def qlucore(
    samples: Samples,
    config: Config,
    label: str,
    logger: LoggerAdapter,
    root: Path,
    workdir: Path,
    executor: Executor,
    **_,
) -> None:
    """Run nf-core/rnaseq (Mapping for qlucore)."""

    if config.qlucore.skip:
        if not config.copy_skipped:
            samples.output = set()
        return samples

    sample_sheet = samples.nfcore_samplesheet(
        location=workdir,
        strandedness=config.strandedness,
    )

    with open(workdir / "nextflow.config", "w") as f:
        if "config" in config.nextflow:
            f.write(f"includeConfig '{config.nextflow.config}'\n\n")
        f.write(qlucore_nf_config)

    (workdir / "dummy.fa").touch()
    (workdir / "dummy.gtf").touch()
    (workdir / "dummy.refflat").touch()

    result, _ = nextflow(
        config.qlucore.nf_main,
        "--starfusion",
        "--skip_qc",
        "--skip_vis",
        "--star_ignore_sjdbgtf",
        f"--fasta {workdir / 'dummy.fa'}",
        f"--transcript {workdir / 'dummy.fa'}",
        f"--gtf {workdir / 'dummy.gtf'}",
        f"--chrgtf {workdir / 'dummy.gtf'}",
        f"--refflat {workdir / 'dummy.refflat'}",
        f"--outdir {workdir}",
        f"--input {sample_sheet}",
        f"--starfusion_ref {config.qlucore.starfusion_ref}",
        f"--starindex_ref {config.qlucore.starfusion_ref}/ref_genome.fa.star.idx",
        f"--read_length {config.read_length}",
        nxf_config=workdir / "nextflow.config",
        config=config,
        name=label,
        workdir=workdir,
        executor=executor,
    )
    result.get()

    for sample in samples:
        executor.submit(
            str(root / "scripts" / "qlucore_subsample.sh"),
            name=f"qlucore_subsample_{sample.id}",
            workdir=workdir,
            cpus=config.qlucore.subsample_threads,
            env={
                "_QLUCORE_SUBSAMPLE_INIT": config.qlucore.subsample_init,
                "_QLUCORE_SUBSAMPLE_FRAC": config.qlucore.subsample_frac,
                "_QLUCORE_SUBSAMPLE_THREADS": config.qlucore.subsample_threads,
                "_QLUCORE_SUBSAMPLE_INPUT_BAM": (
                    workdir
                    / "star_for_starfusion"
                    / f"{sample.id}.Aligned.sortedByCoord.out.bam"
                ),
            },
            callback=partial(
                _subsample_callback,
                logger=logger,
                workdir=workdir,
                sample=sample,
            ),
            error_callback=partial(
                _subsample_error_callback,
                logger=logger,
                sample=sample,
            ),
        )

    executor.wait()

    return samples
