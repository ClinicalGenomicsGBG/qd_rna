"""Module for running nf-core/rnafusion."""

from logging import LoggerAdapter
from pathlib import Path
from copy import deepcopy

from cellophane import cfg, modules
from modules.slims import SlimsSamples
from modules.nextflow import nextflow


def get_output(outdir: Path, config: cfg.Config):
    """Return a dictionary of output files for the rnafusion module."""
    return {
        "arriba": [*(outdir / config.rnaseq.aligner).glob("*.tsv")],
        "arrriba_visualisation": [*(outdir / "arriba_visualisation").glob("*.pdf")],
        "fusioncatcher": [*(outdir / "fusioncatcher").glob("*.txt")],
        "fusioninspector": [*(outdir / "fusioninspector").glob("*")],
        "fusionreport": [*(outdir / "fusionreport").glob("**/*.html")],
        "pizzly": [*(outdir / "pizzly").glob("*")],
        "squid": [*(outdir / "squid").glob("*.txt")],
        "starfusion": [*(outdir / "starfusion").glob("*.tsv")],
    }


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
                f"{sample.id}_{sample.run}" if sample.run else sample.id
            )

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
            name="rnafusion",
            workdir=outdir / "work",
            stderr=outdir / "logs" / f"rnaseq.{timestamp}.err",
            stdout=outdir / "logs" / f"rnaseq.{timestamp}.out",
            cwd=outdir,
        )

        output = get_output(outdir, config)
        logger.debug(output)
