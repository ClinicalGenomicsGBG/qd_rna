"""Base classes and functions for cellophane modules."""

import multiprocessing as mp
import os
import sys
from copy import deepcopy, copy
from dataclasses import dataclass
import logging
from signal import SIGTERM, signal
from typing import Callable, Optional, ClassVar
from pathlib import Path
from queue import Queue

import xxhash
import inspect
import pickle
import psutil

from . import cfg, data, logs


def _cleanup(
    logger: logging.LoggerAdapter,
    output: mp.Queue,
    samples: data.Samples,
) -> Callable:
    def inner(*_):
        for proc in psutil.Process().children(recursive=True):
            logger.debug(f"Waiting for {proc.name()} ({proc.pid})")
            proc.terminate()
            proc.wait()
        output.put(samples)
        output.close()
        output.join_thread()
        raise SystemExit(1)

    return inner


class Runner(mp.Process):
    """Base class for cellophane runners."""

    label: ClassVar[str]
    name: ClassVar[str]
    individual_samples: ClassVar[bool]
    wait: ClassVar[bool]

    def __init_subclass__(
        cls,
        name: str,
        label: Optional[str] = None,
        individual_samples: bool = False,
    ) -> None:
        cls.label = label or name
        cls.name = name
        cls.individual_samples = individual_samples
        super().__init_subclass__()

    def __init__(
        self,
        config: cfg.Config,
        samples: data.Samples,
        timestamp: str,
        log_queue: Queue,
        log_level: int,
        output: mp.Queue,
        root: Path,
    ):
        self.output = output
        self.log_queue = log_queue
        self.log_level = log_level
        super().__init__(
            target=self._main,
            kwargs={
                "config": config,
                "samples": samples,
                "timestamp": timestamp,
                "root": root,
            },
        )

    def _main(
        self,
        config: cfg.Config,
        samples: data.Samples[data.Sample],
        timestamp: str,
        root: Path,
    ) -> Optional[data.Samples[data.Sample]]:
        for sample in samples:
            sample.complete = False
            sample.runner = self.label
        original = deepcopy(samples)

        logger = logs.get_logger(
            label=self.label,
            level=self.log_level,
            queue=self.log_queue,
        )

        signal(SIGTERM, _cleanup(logger, self.output, original))
        sys.stdout = open(os.devnull, "w", encoding="utf-8")
        sys.stderr = open(os.devnull, "w", encoding="utf-8")

        run_hash = xxhash.xxh3_128()
        run_hash.update(pickle.dumps(config))
        run_hash.update(pickle.dumps(inspect.getsource(self.main)))

        for sample in samples:
            for fq in sample.fastq_paths:
                run_hash.update(Path(fq).stat().st_size.to_bytes(8, "big"))
                with open(fq, "rb") as handle:
                    run_hash.update(handle.read(int(128e6)))
                    handle.seek(-int(128e6), 2)
                    run_hash.update(handle.read(int(128e6)))
        
        run_hash = run_hash.hexdigest()[:16]
        outdir: Path = config.outdir / f"{self.name}_{run_hash}"
        logger.debug(f"Run hash: {run_hash}")

        try:
            checksum = pickle.loads((outdir / ".cellophane_integrity").read_bytes())
            integrity = xxhash.xxh3_128()
            integrity.update(pickle.dumps(outdir.glob("**/*")))

            if checksum == integrity.hexdigest():
                returned = pickle.loads((outdir / ".cellophane_cache").read_bytes())
                logger.info("Using cached results")
            else:
                raise Exception("Integrity check failed")

        except Exception as exception:
            logger.debug(f"Unable to load cached results: {exception}")
            returned = None
        
        try:
            returned = returned or self.main(
                samples=samples,
                config=config,
                timestamp=timestamp,
                label=self.label,
                logger=logger,
                root=root,
                outdir=outdir,
            )

            match returned:
                case None:
                    if returned is None and any(
                        s.id not in [o.id for o in original] for s in samples
                    ):
                        logger.warning("Samples were modified but not returned")
                    for sample in original:
                        sample.complete = True
                    self.output.put(original)
                case returned if issubclass(type(returned), data.Samples):
                    for sample in returned:
                        sample.complete = sample.complete or True
                    self.output.put(returned)
                case _:
                    logger.warning(f"Unexpected return type {type(returned)}")
                    self.output.put(original)

            self.output.close()
            self.output.join_thread()

        except Exception as exception:
            logger.critical(exception, exc_info=config.log_level == "DEBUG")
            self.output.put(original)
            self.output.close()
            self.output.join_thread()
            raise SystemExit(1)

    @staticmethod
    def main(**_) -> Optional[data.Samples[data.Sample]]:
        """Main function for the runner."""
        raise NotImplementedError


@dataclass
class Hook:
    """Base class for cellophane pre/post-hooks."""

    label: str
    name: str
    func: Callable
    overwrite: bool
    when: str
    priority: int | float = float("inf")

    def __call__(
        self,
        samples: data.Samples,
        config: cfg.Config,
        timestamp: str,
        log_queue: mp.Queue,
        log_level: int,
        root: Path,
    ) -> data.Samples:
        _logger = logs.get_logger(
            label=self.label,
            level=log_level,
            queue=log_queue,
        )
        return self.func(
            samples=samples,
            config=config,
            timestamp=timestamp,
            logger=_logger,
            root=root,
        )


def pre_hook(
    label: Optional[str] = None,
    overwrite: bool = False,
    priority: int | float = float("inf"),
):
    """Decorator for hooks that will run before all runners."""

    def wrapper(func):
        return Hook(
            label=label or func.__name__,
            name=func.__name__,
            func=func,
            overwrite=overwrite,
            when="pre",
            priority=priority,
        )

    return wrapper


def post_hook(
    label: Optional[str] = None,
    overwrite: bool = False,
    priority: int | float = float("inf"),
):
    """Decorator for hooks that will run after all runners."""

    def wrapper(func):
        return Hook(
            label=label or func.__name__,
            name=func.__name__,
            func=func,
            overwrite=overwrite,
            when="post",
            priority=priority,
        )

    return wrapper


def runner(
    label: Optional[str] = None,
    individual_samples: bool = False,
):
    """Decorator for runners."""

    def wrapper(func):
        class _runner(
            Runner,
            name=func.__name__,
            label=label,
            individual_samples=individual_samples,
        ):
            def __init__(self, *args, **kwargs):
                self.main = staticmethod(func)
                super().__init__(*args, **kwargs)

        return _runner

    return wrapper
