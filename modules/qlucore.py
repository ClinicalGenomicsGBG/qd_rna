"""Module for running nf-core/rnaseq."""

from logging import LoggerAdapter
from pathlib import Path

from cellophane import cfg, modules, data

from modules.qd_rna import Output
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
    withName:NFCORE_RNAFUSION:RNAFUSION:STARFUSION_WORKFLOW:STAR_FOR_STARFUSION {
        ext.args = "--outFilterMultimapNmax 200"
    }
}
"""

@modules.runner(individual_samples=True)
def qlucore(
    samples: data.Samples,
    config: cfg.Config,
    timestamp: str,
    label: str,
    logger: LoggerAdapter,
    root: Path,
    outdir: Path,
) -> data.Samples:
    """Run nf-core/rnaseq (Mapping for qlucore)."""

    if "qlucore" in config and not config.qlucore.skip:
        logger.info("Running STAR + STAR-Fusion for qlucore")

        sample_sheet = samples.nfcore_samplesheet(
            location=outdir,
            strandedness=config.qlucore.strandedness,
        )

        if "workdir" in config.nextflow:
            config.nextflow.workdir.mkdir(parents=True, exist_ok=True)

        with open(outdir / "nextflow.config", "w") as f:
            if "config" in config.nextflow:
                f.write(f"includeConfig {config.qlucore.config}\n\n")
            f.write(qlucore_nf_config)

        nextflow(
            config.qlucore.nf_main,
            "--starfusion",
            "--skip_qc",
            "--skip_vis",
            "--star_ignore_sjdbgtf",
            f"--outdir {outdir}",
            f"--input {sample_sheet}",
            f"--starfusion_ref {config.qlucore.starfusion_ref}",
            f"--read_length {config.qlucore.read_length}",
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

        for sample in samples:
            with open(outdir / f"{sample.id}.qlucore.txt", "w") as f:
                f.write(qlucore_data.format(id=sample.id, run=sample.run or ""))

            sample.output = [
                Output(
                    src = outdir.glob(f"{sample.id}.qlucore.txt"),
                    dest_dir = Path(sample.id) / "qlucore",
                ),
                Output(
                    src = (outdir / "star_for_starfusion").glob(f"{sample.id}.*.bam"),
                    dest_dir = Path(sample.id) / "qlucore",
                ),
                Output(
                    src = (outdir / "starfusion").glob(f"{sample.id}.*.tsv"),
                    dest_dir = Path(sample.id) / "qlucore",
                ),
            ]
        
    return samples
