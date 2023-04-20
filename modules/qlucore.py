"""Module for running nf-core/rnaseq."""

from multiprocessing import Process
from logging import LoggerAdapter
from pathlib import Path

from cellophane import cfg, modules, data, sge

from modules.nextflow import nextflow


qlucore_data = """\
<Header Producer='Qlucore' Format='PatientData' FormatVersion='0.1' QFFVersion='1.0'/>
<PatientData>
  <PatientName>PATIENT NAME</PatientName>
  <PatientId>{id}</PatientId>
  <SampleOrigin>{run}</SampleOrigin>
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


@modules.runner()
def qlucore(
    samples: data.Samples,
    config: cfg.Config,
    timestamp: str,
    label: str,
    logger: LoggerAdapter,
    root: Path,
    outdir: Path,
) -> None:
    """Run nf-core/rnaseq (Mapping for qlucore)."""

    if not config.qlucore.skip:
        sample_sheet = samples.nfcore_samplesheet(
            location=outdir,
            strandedness=config.strandedness,
        )

        if "workdir" in config.nextflow:
            config.nextflow.workdir.mkdir(parents=True, exist_ok=True)

        with open(outdir / "nextflow.config", "w") as f:
            if "config" in config.nextflow:
                f.write(f"includeConfig '{config.nextflow.config}'\n\n")
            f.write(qlucore_nf_config)

        (outdir / "dummy.fa").touch()
        (outdir / "dummy.gtf").touch()
        (outdir / "dummy.refflat").touch()

        nextflow(
            config.qlucore.nf_main,
            "--starfusion",
            "--skip_qc",
            "--skip_vis",
            "--star_ignore_sjdbgtf",
            f"--fasta {outdir / 'dummy.fa'}",
            f"--transcript {outdir / 'dummy.fa'}",
            f"--gtf {outdir / 'dummy.gtf'}",
            f"--chrgtf {outdir / 'dummy.gtf'}",
            f"--refflat {outdir / 'dummy.refflat'}",
            f"--outdir {outdir}",
            f"--input {sample_sheet}",
            f"--starfusion_ref {config.qlucore.starfusion_ref}",
            f"--starindex_ref {config.qlucore.starfusion_ref}/ref_genome.fa.star.idx",
            f"--read_length {config.read_length}",
            config=config,
            name=label,
            nf_config=outdir / "nextflow.config",
            workdir=outdir / "work",
            report=outdir / "logs" / f"{label}.{timestamp}.nextflow_report.html",
            log=outdir / "logs" / f"{label}.{timestamp}.nextflow.log",
            stderr=outdir / "logs" / f"{label}.{timestamp}.nextflow.err",
            stdout=outdir / "logs" / f"{label}.{timestamp}.nextflow.out",
            cwd=outdir,
        )

        _SUBSAMPLE_PROCS: dict[str, Process] = {}
        for sample in samples:
            _SUBSAMPLE_PROCS[sample.id] = sge.submit(
                str(root / "scripts" / "qlucore_subsample.sh"),
                queue=config.nextflow.sge_queue,
                name=f"qlucore_subsample_{sample.id}",
                pe=config.nextflow.sge_pe,
                slots=config.qlucore.subsample_threads,
                stdout=outdir / "logs" / f"{sample.id}.qlucore_subsample.out",
                stderr=outdir / "logs" / f"{sample.id}.qlucore_subsample.err",
                cwd=outdir,
                check=False,
                env={
                    "_QLUCORE_SUBSAMPLE_FRAC": config.qlucore.subsample_frac,
                    "_QLUCORE_SUBSAMPLE_THREADS": config.qlucore.subsample_threads,
                    "_QLUCORE_SUBSAMPLE_INPUT_BAM": (
                        outdir
                        / "star_for_starfusion"
                        / f"{sample.id}.Aligned.sortedByCoord.out.bam"
                    ),
                },
            )

        for sample in samples:
            subsample_proc = _SUBSAMPLE_PROCS[sample.id]
            subsample_proc.join()
            if subsample_proc.exitcode != 0:
                logger.error(f"Subsampling failed for {sample.id}")

            with open(outdir / f"{sample.id}.qlucore.txt", "w") as f:
                f.write(qlucore_data.format(id=sample.id, run=sample.run or ""))

    if not config.qlucore.skip or config.rsync.copy_skipped:
        for sample in samples:
            sample.output = [
                data.Output(
                    src=outdir.glob(f"{sample.id}.qlucore.txt"),
                    dest_dir=Path(sample.id) / "qlucore",
                ),
                data.Output(
                    src=(outdir / "star_for_starfusion").glob(f"{sample.id}.*.bam"),
                    dest_dir=Path(sample.id) / "qlucore",
                ),
                data.Output(
                    src=(outdir / "starfusion").glob(f"{sample.id}.*.tsv"),
                    dest_dir=Path(sample.id) / "qlucore",
                ),
            ]

    return samples
