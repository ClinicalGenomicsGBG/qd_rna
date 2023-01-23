"""Base classes and functions for cellophane modules."""

import multiprocessing as mp
import os
import sys
from copy import deepcopy
from dataclasses import dataclass
import logging
from signal import SIGTERM, signal
from typing import Callable, Optional, ClassVar
from pathlib import Path
from queue import Queue
from functools import cached_property

import psutil

from . import cfg, data, logs


def _cleanup(logger: logging.LoggerAdapter):
    def inner(*_):
        for proc in psutil.Process().children(recursive=True):
            logger.debug(f"Waiting for {proc.name()} ({proc.pid})")
            proc.terminate()
            proc.wait()
        raise SystemExit(1)

    return inner


class Runner(mp.Process):
    """Base class for cellophane runners."""

    label: ClassVar[str]
    individual_samples: ClassVar[bool]
    wait: ClassVar[bool]

    def __init_subclass__(
        cls,
        label: Optional[str] = None,
        individual_samples: bool = False,
    ) -> None:
        cls.label = label or cls.__name__
        cls.individual_samples = individual_samples
        super().__init_subclass__()

    def __init__(
        self,
        config: cfg.Config,
        samples: data.Samples,
        log_queue: Queue,
        log_level: int,
        output_queue: Queue,
        root: Path,
    ):
        super().__init__(
            target=self._main,
            kwargs={
                "label": self.label,
                "config": config,
                "samples": samples,
                "log_queue": log_queue,
                "log_level": log_level,
                "output_queue": output_queue,
                "root": root,
            },
        )

    def _main(
        self,
        label: str,
        config: cfg.Config,
        samples: data.Samples,
        log_queue: Queue,
        log_level: int,
        output_queue: Queue,
        root: Path,
    ) -> None:
        logger = logs.get_logger(
            label=label,
            level=log_level,
            queue=log_queue,
        )
        signal(SIGTERM, _cleanup(logger))
        sys.stdout = open(os.devnull, "w", encoding="utf-8")
        sys.stderr = open(os.devnull, "w", encoding="utf-8")
        try:
            original = deepcopy(samples)
            returned = self.main(
                samples=samples,
                config=config,
                label=label,
                logger=logger,
                root=root,
            )

            match returned:
                case None if any(s not in original for s in samples):
                    logger.warning("Runner returned None, but samples were modified")
                    output_queue.put(original)
                case data.Samples | None:
                    output_queue.put(returned)
                case _:
                    logger.warning(
                        f"Runner returned an unexpected type {type(returned)}"
                    )
                    output_queue.put(original)

        except Exception as exception:
            logger.critical(
                "Caught an exception",
                exc_info=config.log_level <= logging.DEBUG,
            )
            raise SystemExit(1) from exception

    @cached_property
    def processed_samples(self) -> data.Samples | None:
        """Return the processed samples of the runner."""
        return self._output_queue.get_nowait()

    @staticmethod
    def main(*args, **kwargs) -> None:
        """Main function for the runner."""
        raise NotImplementedError


@dataclass
class Hook:
    """Base class for cellophane pre/post-hooks."""

    label: str
    func: Callable
    overwrite: bool
    when: str
    priority: int | float = float("inf")

    def __call__(
        self,
        config: cfg.Config,
        samples: data.Samples,
        log_queue: mp.Queue,
        log_level: int,
        root: Path,
    ) -> data.Samples:
        _adapter = logs.get_logger(
            label=self.label,
            level=log_level,
            queue=log_queue,
        )
        return self.func(
            config=config,
            samples=samples,
            logger=_adapter,
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
            label=label or func.__name__,
            individual_samples=individual_samples,
        ):
            @staticmethod
            def main(*args, **kwargs):
                return func(*args, **kwargs)

        return _runner

    return wrapper
