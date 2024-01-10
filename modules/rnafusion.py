"""Module for running nf-core/rnafusion."""

from pathlib import Path

from cellophane import output, runner, Executor, Config, Sample, Samples
from modules.nextflow import nextflow

@output(
    "arriba_visualisation/{sample.id}.pdf",
    dst_dir="{sample.id}/arriba",
)
@output(
    "arriba/{sample.id}.*.tsv",
    dst_dir="{sample.id}/arriba",
)
@output(
    "fusioncatcher/{sample.id}.*.txt",
    dst_dir="{sample.id}/fusioncatcher",
)
@output(
    "fusioninspector/{sample.id}.*",
    dst_dir="{sample.id}/fusioninspector",
)
@output(
    "fusionreport/{sample.id}",
    dst_dir="{sample.id}/fusionreport",
)
@output(
    "pizzly/{sample.id}.*",
    dst_dir="{sample.id}/pizzly",
)
@output(
    "squid/{sample.id}.*.txt",
    dst_dir="{sample.id}/squid",
)
@output(
    "star_for_starfusion/{sample.id}.*.ba*",
    dst_dir="{sample.id}",
)
@runner()
def rnafusion(
    samples: Samples,
    config: Config,
    label: str,
    workdir: Path,
    executor: Executor,
    **_,
) -> Samples:
    """Run nf-core/rnafusion."""
    if config.rnafusion.skip:
        if not config.copy_skipped:
            samples.output = set()
        return samples
    
    sample_sheet = samples.nfcore_samplesheet(
        location=workdir,
        strandedness=config.strandedness,
    )

    result, _ = nextflow(
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

    result.get()

    return samples
