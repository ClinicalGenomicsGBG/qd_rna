"""Module for fetching files from HCP."""

from pathlib import Path
from typing import Optional, Mapping

from cellophane import cfg, sge, modules, data


class NextflowSamples(data.Mixin):
    def nfcore_samplesheet(self, *_, location: str | Path, **kwargs) -> Path:
        """Write a Nextflow samplesheet"""
        Path(location).mkdir(parents=True, exist_ok=True)
        _data = [
            {
                "sample": sample.id,
                "fastq_1": str(sample.fastq_paths[0]),
                "fastq_2": str(sample.fastq_paths[1]),
                **{
                    k: v[sample.id] if isinstance(v, Mapping) else v
                    for k, v in kwargs.items()
                },
            }
            for sample in self
        ]

        _header = ",".join(_data[0].keys())

        _samplesheet = "\n".join([_header, *(",".join(d.values()) for d in _data)])
        _path = Path(location) / "samples.nextflow.csv"
        with open(_path, "w", encoding="utf-8") as handle:
            handle.write(_samplesheet)

        return _path

def nextflow(
    main: Path,
    *args,
    config: cfg.Config,
    env: dict[str, str] = {},
    log: Optional[Path] = None,
    report: Optional[Path] = None,
    workdir: Optional[Path] = None,
    nf_config: Optional[Path] = None,
    ansi_log: bool = False,
    resume: bool = False,
    **kwargs,
):
    if "workdir" in config.nextflow:
        config.nextflow.workdir.mkdir(parents=True, exist_ok=True)

    sge.submit(
        str(Path(__file__).parent / "scripts" / "nextflow.sh"),
        f"-log {log}" if log else "",
        (
            f"-config {nf_config}"
            if nf_config
            else f"-config {config.nextflow.config}"
            if "config" in config.nextflow
            else ""
        ),
        f"run {main}",
        "-ansi-log false" if not ansi_log or config.nextflow.ansi_log else "",
        f"-work-dir {workdir}" if workdir else "",
        "-resume" if resume else "",
        f"-with-report {report}" if report else "",
        f"-profile {config.nextflow.profile}",
        *args,
        env={
            "_NXF_MODULE": config.nextflow.nf_module,
            "_JAVA_MODULE": config.nextflow.java_module,
            **env,
        },
        queue=config.nextflow.sge_queue,
        pe=config.nextflow.sge_pe,
        slots=config.nextflow.sge_slots,
        **kwargs,
    )
