from cellophane import Checkpoints, Config, Executor, Samples, Sample, output, runner
from pathlib import Path
from logging import LoggerAdapter
from xml.etree import ElementTree as ET
from datetime import datetime
from functools import partial


def generate_qsd(sample: Sample, outpath: Path) -> None:
    """Generate a Qlucore Sample Data (QSD) XML file for the given sample."""
    qff = ET.Element(
        "QFF",
        {
            "Producer": "Qlucore",
            "Format": "QlucoreSampleData",
            "FormatVersion": "1.0",
            "QFFVersion": "1.1",
        },
    )

    sample_data = ET.SubElement(qff, "SampleData")
    ET.SubElement(sample_data, "SubjectId").text = sample.id
    ET.SubElement(sample_data, "SubjectName").text = sample.id
    ET.SubElement(sample_data, "SampleDateTime").text = datetime.now().strftime(
        "%Y-%b-%d %H:%M:%S"
    )
    ET.SubElement(sample_data, "SampleId").text = sample.id
    ET.SubElement(sample_data, "SampleTissue").text = "Blood sample"

    xml_content = ET.tostring(qff, encoding="unicode", xml_declaration=True)
    with open(outpath, "w") as f:
        f.write(xml_content)


def _calculate_samtools_fraction(star_log: Path, target_reads: int, logger: LoggerAdapter) -> float | None:
    """Calculate the fraction for samtools subsampling based on STAR Log.final.out."""
    if not star_log.exists():
        logger.error(f"STAR log file not found: {star_log}")
        return None
    log_text = star_log.read_text()
    match = re.search(r"Uniquely mapped reads number[^\d]+(\d+)", log_text)
    if not match:
        logger.error(f"Could not find uniquely mapped reads in STAR log: {star_log}")
        return None
    uniquely_mapped = int(match.group(1))
    if uniquely_mapped <= 0:
        logger.error(f"No uniquely mapped reads found in STAR log: {star_log}")
        return None
    fraction = target_reads / uniquely_mapped
    return min(fraction, 1.0)  # Cap at 1.0

def _star_error_callback(
    samples: Samples, exc: Exception, logger: LoggerAdapter
) -> None:
    reason = f"STAR failed for {samples[0].id}: {exc}"
    logger.error(reason)
    for sample in samples:
        sample.fail(reason)


def _samtools_error_callback(
    samples: Samples, exc: Exception, logger: LoggerAdapter
) -> None:
    reason = f"Samtools subsampling failed for {samples[0].id}: {exc}"
    logger.error(reason)
    for sample in samples:
        sample.fail(reason)


DEFAULT_STAR_ARGS = [
    "--outSAMattrRGline", "ID:GRPundef",
    "--twopassMode", "Basic",
    "--outReadsUnmapped", "None",
    "--readFilesCommand", "zcat",
    "--outSAMtype", "BAM", "SortedByCoordinate",
    "--outSAMstrandField", "intronMotif",
    "--outSAMunmapped", "Within",
    "--chimSegmentMin", "12",
    "--chimJunctionOverhangMin", "8",
    "--chimOutJunctionFormat", "1",
    "--alignSJDBoverhangMin", "10",
    "--alignMatesGapMax", "100000",
    "--alignIntronMax", "100000",
    "--alignSJstitchMismatchNmax", "5", "-1", "5", "5",
    "--chimMultimapScoreRange", "3",
    "--chimScoreJunctionNonGTAG", "-4",
    "--chimMultimapNmax", "20",
    "--chimNonchimScoreDropMin", "10",
    "--peOverlapNbasesMin", "12",
    "--peOverlapMMp", "0.1",
    "--alignInsertionFlush", "Right",
    "--alignSplicedMateMapLminOverLmate", "0",
    "--alignSplicedMateMapLmin", "30",
    "--outFilterMultimapNmax", "200",
]


@output(
    "Aligned.sortedByCoord.out.bam",
    dst_dir="{sample.id}_{sample.last_run}_%y%m%d-%H%M%S/qlucore",
    checkpoint="star",
)
@output(
    "subsampled.bam",
    dst_dir="{sample.id}_{sample.last_run}_%y%m%d-%H%MS/qlucore",
    checkpoint="subsample",
)
@output(
    "{sample.id}.qsd",
    dst_dir="{sample.id}_{sample.last_run}_%y%m%d-%H%M%S/qlucore"
)
@runner(split_by="id")
def qlucore(
    samples: Samples,
    config: Config,
    logger: LoggerAdapter,
    root: Path,
    workdir: Path,
    executor: Executor,
    checkpoints: Checkpoints,
    **_,
) -> None:
    """Run STAR + samtools subsampling for qlucore."""
    if not checkpoints.star.check():
        fw_reads = [sample.files[0] for sample in samples]
        rw_reads = [sample.files[1] for sample in samples]
        bind_paths = set(
            [str(Path(config.qlucore.star.index).parent)]
            + [str(file.parent) for file in fw_reads + rw_reads]
        )
        bind_args = ["--bind", ",".join(bind_paths)] if bind_paths else []
        star_result, star_uuid = executor.submit(
            "apptainer exec",
            *bind_args,
            config.qlucore.star.container,
            "STAR",
            "--readFilesIn", ",".join(fw_reads), ",".join(rw_reads),
            "--runThreadN", config.qlucore.star.threads,
            "--genomeDir", config.qlucore.star.index,
            *DEFAULT_STAR_ARGS,
            workdir=workdir,
            cpus=config.qlucore.star.threads,
            name=f"STAR_{samples[0].id}",
            wait=True,
            error_callback=partial(
                _star_error_callback, samples=samples, logger=logger
            ),
        )
        executor.wait(star_uuid)

    if not checkpoints.subsample.check():
        subsample_fraction = _calculate_samtools_fraction(
            workdir / "Log.final.out",
            target_reads=config.qlucore.samtools.target,
            logger=logger
        ) or config.qlucore.samtools.fraction
        if subsample_fraction == 1.0:
            logger.info("No subsampling needed for %s (fraction=1.0)", samples[0].id)
            # Just Crete a symlink to avoid unnecessary work
            (workdir / "subsampled.bam").symlink_to(workdir / "Aligned.sortedByCoord.out.bam")
        else:
            samtools_extra_args = config.qclucore.samtools.args or []
            samtools_result, samtools_uuid = executor.submit(
                "appatainer exec",
                config.samtools.container,
                "samtools view",
                "-s", f"{subsample_fraction:.6f}",
                "-@", config.qlucore.samtools.threads,
                "-o", workdir / "subsampled.bam",
                *samtools_extra_args,
                "Aligned.sortedByCoord.out.bam",
                cpus=config.qlucore.samtools.threads,
                name=f"samtools_subsample_{samples[0].id}",
                wait=True,
                error_callback=partial(
                    _samtools_error_callback, samples=samples, logger=logger
                ),
            )
        executor.wait(samtools_uuid)

    if not checkpoints.main.check():
        generate_qsd(samples[0], workdir / f"{samples[0].id}.qsd")
