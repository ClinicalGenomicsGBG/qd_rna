"""Module for running nf-core/rnafusion."""

from logging import LoggerAdapter
from pathlib import Path

from cellophane import cfg, modules, data
from modules.nextflow import nextflow
from modules.qd_rna import Output


@modules.runner()
def rnafusion(
    samples: data.Samples,
    config: cfg.Config,
    timestamp: str,
    label: str,
    logger: LoggerAdapter,
    root: Path,
    outdir: Path,
) -> data.Samples:
    """Run nf-core/rnafusion."""

    if "rnafusion" in config and not config.rnafusion.skip:
        logger.info("Running nf-core/rnafusion")

        sample_sheet = samples.nfcore_samplesheet(
            location=outdir,
            strandedness=config.rnafusion.strandedness,
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
            f"--read_length {config.rnafusion.read_length}",
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

        for sample in samples:
            sample.output = [
                Output(
                    src = [
                        *(outdir / "arriba_visualisation").glob(f"{sample.id}.pdf"),
                        *(outdir / "arriba").glob(f"{sample.id}.*.tsv"),
                    ],
                    dest_dir = Path(sample.id) / "arriba",
                ),
                Output(
                    src = (outdir / "fusioncatcher").glob(f"{sample.id}.*.txt"),
                    dest_dir = Path(sample.id) / "fusioncatcher",
                ),
                Output(
                    src = (outdir / "fusioninspector").glob(f"{sample.id}.*"),
                    dest_dir = Path(sample.id) / "fusioninspector",
                ),
                Output(
                    src = (outdir / "fusionreport").glob("sample.id/*.html"),
                    dest_dir = Path(sample.id) / "fusionreport",
                ),
                Output(
                    src = (outdir / "pizzly").glob(f"{sample.id}.*"),
                    dest_dir = Path(sample.id) / "pizzly",
                ),
                Output(
                    src = (outdir / "squid").glob(f"{sample.id}.*.txt"),
                    dest_dir = Path(sample.id) / "squid",
                ),
                Output(
                    src = (outdir / "starfusion").glob(f"{sample.id}.*.tsv"),
                    dest_dir = Path(sample.id) / "starfusion",
                ),
            ]

    return samples