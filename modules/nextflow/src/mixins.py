"""Nextflow-specific mixins for Cellophane."""

from pathlib import Path

from cellophane import Samples


class NextflowSamples(Samples):
    """Samples with Nextflow-specific methods."""

    def nfcore_samplesheet(
        self,
        *_,
        location: str | Path,
        **kwargs,
    ) -> Path:
        """Write a Nextflow samplesheet"""
        Path(location).mkdir(parents=True, exist_ok=True)
        _data = [
            {
                "sample": sample.id,
                "fastq_1": str(sample.files[0]),
                "fastq_2": str(sample.files[1]) if len(sample.files) > 1 else "",
                **{
                    k: (
                        v.format(sample=sample)
                        if isinstance(v, str)
                        else v
                    )
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
