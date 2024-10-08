from logging import LoggerAdapter
from pathlib import Path
from typing import Callable, Iterator

from cellophane import cfg, executors
from mpire.async_result import AsyncResult


class Extractor:
    """Base class for extractors."""

    label: str
    script: Path
    suffixes: tuple[str, ...]
    conda_spec: dict | None

    def __init_subclass__(
        cls,
        label: str,
        script: Path,
        suffixes: tuple[str, ...],
        conda_spec: dict | None = None,
    ) -> None:
        cls.label = label
        cls.script = script
        cls.suffixes = suffixes
        cls.conda_spec = conda_spec

    def basename(self, path: Path) -> str:
        _name = path.name
        for suffix in self.suffixes:
            _name = _name.replace(suffix, "")
        return _name

    def extracted_paths(
        self,
        workdir: Path,
        compressed_path: Path,
    ) -> Iterator[Path]:  # pragma: no cover
        del workdir, compressed_path  # Unused
        raise NotImplementedError

    def extract(
        self,
        *,
        logger: LoggerAdapter,
        compressed_path: Path,
        workdir: Path,
        config: cfg.Config,
        env: dict | None = None,
        executor: executors.Executor,
        callback: Callable | None = None,
        error_callback: Callable | None = None,
    ) -> AsyncResult | None:
        """
        Summary:
        Extracts a compressed file using the provided executor and configuration.

        Args:
            logger: LoggerAdapter - The logger to use for logging messages.
            compressed_path: Path - The path to the compressed file to extract.
            config: cfg.Config - The configuration settings to use for extraction.
            env: dict | None - Optional environment variables to set during extraction.
            executor: executors.Executor - The executor to use for jobs.
            callback: Callable | None - Optional callback called on completion.
            error_callback: Callable | None - Optional callback called on error.

        Returns:
            AsyncResult | None - The result of the extraction process,
                or None if already extracted.
        """
        if not compressed_path.exists():
            logger.error(f"Compressed file {compressed_path.name} not found")
            return None

        if [*self.extracted_paths(workdir, compressed_path)]:  # pylint: disable=using-constant-test
            logger.debug(f"Already extracted {compressed_path.name}")
            if callback is not None:
                callback(None)
            return None

        else:
            logger.info(f"Extracting {compressed_path.name} with {self.label}")
            result, _ = executor.submit(
                self.script,
                name=f"unpack_{compressed_path.name}",
                env={
                    **(env or {}),
                    "UNPACK_INIT": config.unpack.get("init", ""),
                    "UNPACK_EXIT": config.unpack.get("exit", ""),
                    "COMPRESSED_PATH": compressed_path,
                    "EXTRACTED_PATH": workdir / f"{self.basename(compressed_path)}.fastq.gz",
                    "THREADS": config.unpack.threads,
                },
                cpus=config.unpack.threads,
                workdir=workdir,
                callback=callback,
                error_callback=error_callback,
                conda_spec=self.conda_spec,
            )

            return result


class PetageneExtractor(
    Extractor,
    label="petagene",
    script=Path(__file__).parent.parent / "scripts" / "petagene.sh",
    suffixes=(".fasterq",),
):
    """Petagene extractor."""

    def extracted_paths(self, workdir: Path, compressed_path: Path) -> Iterator[Path]:
        _base = self.basename(compressed_path)
        if (extracted := workdir / f"{_base}.fastq.gz").exists():
            yield extracted


class SpringExtractor(
    Extractor,
    label="spring",
    script=Path(__file__).parent.parent / "scripts" / "spring.sh",
    suffixes=(".spring",),
    conda_spec={"dependencies": ["spring >=1.1.1, <2.0.0"]},
):
    """Spring extractor."""

    def extracted_paths(self, workdir: Path, compressed_path: Path) -> Iterator[Path]:
        _base = self.basename(compressed_path)
        if (extracted := workdir / f"{_base}.fastq.gz").exists():
            yield extracted

        elif all(
            (
                (extracted1 := workdir / f"{_base}.1.fastq.gz").exists(),
                (extracted2 := workdir / f"{_base}.2.fastq.gz").exists(),
            )
        ):
            yield extracted1
            yield extracted2

        elif all(
            (
                (extracted1 := workdir / f"{_base}.fastq.gz.1").exists(),
                (extracted2 := workdir / f"{_base}.fastq.gz.2").exists(),
            )
        ):
            yield extracted1.rename(workdir / f"{_base}.1.fastq.gz")
            yield extracted2.rename(workdir / f"{_base}.2.fastq.gz")
