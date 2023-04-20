"""Module for running nf-core/rnafusion."""

from pathlib import Path

from cellophane import cfg, modules, data
from modules.nextflow import nextflow


@modules.runner()
def rnafusion(
    samples: data.Samples,
    config: cfg.Config,
    timestamp: str,
    label: str,
    outdir: Path,
    **_,
) -> data.Samples:
    """Run nf-core/rnafusion."""

    if not config.rnafusion.skip:
        sample_sheet = samples.nfcore_samplesheet(
            location=outdir,
            strandedness=config.strandedness,
        )

        if "workdir" in config.nextflow:
            config.nextflow.workdir.mkdir(parents=True, exist_ok=True)

        nextflow(
            config.rnafusion.nf_main,
            f"--outdir {outdir}",
            f"--input {sample_sheet}",
            f"--genomes_base {config.rnafusion.genomes_base}",
            f"--arriba_ref {config.rnafusion.arriba_ref}",
            f"--arriba_ref_blacklist {config.rnafusion.arriba_blacklist}",
            f"--arriba_ref_protein_domain {config.rnafusion.arriba_protein_domain}",
            f"--fusionreport_tool_cutoff {config.rnafusion.fusionreport_tool_cutoff}",
            f"--read_length {config.read_length}",
            "--fusioncatcher_limitSjdbInsertNsj 4000000",
            "--all",
            "--fusioninspector_filter",
            config=config,
            name=label,
            workdir=outdir / "work",
            report=outdir / "logs" / f"{label}.{timestamp}.nextflow_report.html",
            log=outdir / "logs" / f"{label}.{timestamp}.nextflow.log",
            stderr=outdir / "logs" / f"{label}.{timestamp}.nextflow.err",
            stdout=outdir / "logs" / f"{label}.{timestamp}.nextflow.out",
            cwd=outdir,
        )

    if not config.rnafusion.skip or config.results.copy_skipped:
        for sample in samples:
            sample.output = [
                data.Output(
                    src=[
                        *(outdir / "arriba_visualisation").glob(f"{sample.id}.pdf"),
                        *(outdir / "arriba").glob(f"{sample.id}.*.tsv"),
                    ],
                    dest_dir=Path(sample.id) / "arriba",
                ),
                data.Output(
                    src=(outdir / "fusioncatcher").glob(f"{sample.id}.*.txt"),
                    dest_dir=Path(sample.id) / "fusioncatcher",
                ),
                data.Output(
                    src=(outdir / "fusioninspector").glob(f"{sample.id}.*"),
                    dest_dir=Path(sample.id) / "fusioninspector",
                ),
                data.Output(
                    src=outdir / "fusionreport" / sample.id,
                    dest_dir=Path(sample.id) / "fusionreport",
                ),
                data.Output(
                    src=(outdir / "pizzly").glob(f"{sample.id}.*"),
                    dest_dir=Path(sample.id) / "pizzly",
                ),
                data.Output(
                    src=(outdir / "squid").glob(f"{sample.id}.*.txt"),
                    dest_dir=Path(sample.id) / "squid",
                ),
                data.Output(
                    src=(outdir / "star_for_starfusion").glob(f"{sample.id}.*.bam"),
                    dest_dir=Path(sample.id),
                ),
                data.Output(
                    src=(outdir / "samtools_index_for_qc").glob(f"{sample.id}.*.bai"),
                    dest_dir=Path(sample.id),
                ),
                data.Output(
                    src=(outdir / "starfusion").glob(f"{sample.id}.*.tsv"),
                    dest_dir=Path(sample.id) / "starfusion",
                ),
            ]

    return samples
