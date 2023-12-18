"""Module for running nf-core/rnafusion."""

from pathlib import Path

from cellophane import cfg, modules, data, executors
from modules.nextflow import nextflow

@data.output(
    "arriba_visualisation/{sample.id}.pdf",
    dst_dir="{sample.id}/arriba",
)
@data.output(
    "arriba/{sample.id}.*.tsv",
    dst_dir="{sample.id}/arriba",
)
@data.output(
    "fusioncatcher/{sample.id}.*.txt",
    dst_dir="{sample.id}/fusioncatcher",
)
@data.output(
    "fusioninspector/{sample.id}.*",
    dst_dir="{sample.id}/fusioninspector",
)
@data.output(
    "fusionreport/{sample.id}",
    dst_dir="{sample.id}/fusionreport",
)
@data.output(
    "pizzly/{sample.id}.*",
    dst_dir="{sample.id}/pizzly",
)
@data.output(
    "squid/{sample.id}.*.txt",
    dst_dir="{sample.id}/squid",
)
@data.output(
    "star_for_starfusion/{sample.id}.*.ba*",
    dst_dir="{sample.id}",
)
@modules.runner()
def rnafusion(
    samples: data.Samples,
    config: cfg.Config,
    label: str,
    workdir: Path,
    executor: executors.Executor,
    **_,
) -> data.Samples:
    """Run nf-core/rnafusion."""
    if config.rnafusion.skip:
        if not config.copy_skipped:
            samples.output = set()
        return samples
    
    sample_sheet = samples.nfcore_samplesheet(
        location=workdir,
        strandedness=config.strandedness,
    )

    nextflow(
        config.rnafusion.nf_main,
        f"--outdir {workdir}",
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
        workdir=workdir,
        executor=executor,
    )

    return samples
