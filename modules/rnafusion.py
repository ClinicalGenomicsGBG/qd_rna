"""Module for running nf-core/rnafusion."""

from logging import LoggerAdapter
from pathlib import Path
from copy import deepcopy

from cellophane import cfg, modules
from modules.slims import SlimsSamples
from modules.nextflow import nextflow
from modules.qd_rna import Output


@modules.runner()
def rnafusion(
    samples: SlimsSamples,
    config: cfg.Config,
    timestamp: str,
    label: str,
    logger: LoggerAdapter,
    root: Path,
    outdir: Path,
) -> None:
    """Run nf-core/rnafusion."""

    if "rnafusion" in config and not config.rnafusion.skip:
        logger.info("Running nf-core/rnafusion")
        logger.debug(f"Output will be written to {outdir}")

        _samples = deepcopy(samples)
        for idx, sample in enumerate(_samples):
            _samples[idx].id = (
                f"{sample.id}_{sample.run}"
                if "run" in sample and sample.run
                else sample.id
            )

        sample_sheet = _samples.nfcore_samplesheet(
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

        for sample in _samples:
            samples.output += [
                Output(
                    src = [
                        *(outdir / "arriba_visualisation" / f"{sample.id}.pdf"),
                        *(outdir / "arriba" / f"{sample.id}.*.tsv"),
                    ],
                    dest = (Path(sample.id) / "arriba"),
                ),
                Output(
                    src = (outdir / "fusioncatcher").glob(f"{sample.id}.*.txt"),
                    dest = (Path(sample.id) / "fusioncatcher"),
                ),
                Output(
                    src = (outdir / "fusioninspector").glob(f"{sample.id}.*"),
                    dest = (Path(sample.id) / "fusioninspector"),
                ),
                Output(
                    src = (outdir / "fusionreport" / sample.id).glob("*.html"),
                    dest = (Path(sample.id) / "fusionreport"),
                ),
                Output(
                    src = (outdir / "pizzly").glob(f"{sample.id}.*"),
                    dest = (Path(sample.id) / "pizzly"),
                ),
                Output(
                    src = (outdir / "squid").glob(f"{sample.id}.*.txt"),
                    dest = (Path(sample.id) / "squid"),
                ),
                Output(
                    src = (outdir / "starfusion").glob(f"{sample.id}.*.tsv"),
                    dest = (Path(sample.id) / "starfusion"),
                ),
            ]

        return samples
