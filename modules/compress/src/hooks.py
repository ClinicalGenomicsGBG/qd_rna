from logging import LoggerAdapter
from pathlib import Path

from cellophane import Config, Samples, Output, Executor, post_hook

from .util import cram_compress, should_convert, find_bai_output


@post_hook(label="Compress Output", per="sample", before="rsync_results", after="slims_sync_post", condition="complete")
def compress_results(
    samples: Samples,
    logger: LoggerAdapter,
    config: Config,
    workdir: Path,
    executor: Executor,
    **_,
):
    if not getattr(samples, "output", None):
        logger.warning("No output to compress")
        return samples

    if not getattr(config, "compress", None) or not config.compress.enabled:
        logger.info("Compression is disabled")
        return samples

    rules = list(getattr(config.compress, "rules", []) or [])
    threads = int(getattr(config.compress, "threads", 4))

    if not rules:
        logger.info("No files specified for compression")
        return samples

    new_outputs = set(samples.output)  # set of original outputs, will be modified in place

    for out in list(samples.output):
        convert, method, ref = should_convert(out, rules, workdir, logger)
        if not convert:
            continue

        logger.info(f"Compressing {out.src} using method '{method}'")
        if method == "cram":
            bam_src = Path(out.src)
            bam_dst = Path(out.dst)

            cram_src = bam_src.with_suffix(".cram")
            cram_dst = bam_dst.with_suffix(".cram")

            # submit cram job
            cram_compress(
                executor=executor,
                input_bam=bam_src,
                output_cram=cram_src,
                reference=ref,
                threads=threads,
                logger=logger,
                name=f"cram_{bam_src.stem}",
                workdir=workdir / "compress" / bam_src.stem,
            )

            # swap BAM -> CRAM
            new_outputs.discard(out)
            new_outputs.add(Output(
                src=cram_src,
                dst=cram_dst,
                optional=out.optional,
                checkpoint=out.checkpoint,
            ))

            # if there's an Output for BAM index, replace it with CRAI
            bai_out = find_bai_output(new_outputs, bam_src)
            if bai_out is not None:
                new_outputs.discard(bai_out)

                crai_src = cram_src.with_suffix(".cram.crai")
                crai_dst = cram_dst.with_suffix(".cram.crai")

                new_outputs.add(Output(
                    src=crai_src,
                    dst=crai_dst,
                    optional=True,
                    checkpoint=bai_out.checkpoint, 
                ))

        else:
            logger.warning(f"Unsupported compression method '{method}' for {out.src}, skipping")

    executor.wait()  # wait for all compression jobs to finish before proceeding

    samples.output = new_outputs

    logger.debug("COMPRESS hook: samples id=%s output id=%s", id(samples), id(samples.output))
    logger.debug("COMPRESS hook: any cram=%s any bam=%s",
                any(str(o.src).endswith(".cram") for o in samples.output),
                any(str(o.src).endswith(".bam") for o in samples.output))

    return samples

