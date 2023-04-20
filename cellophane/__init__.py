"""Cellophane: A library for writing modular wrappers"""
import inspect
import logging
import multiprocessing as mp
import time
from copy import deepcopy
from pathlib import Path
from typing import Any, Callable, Optional
from uuid import UUID

import rich_click as click
import yaml
from humanfriendly import format_timespan
from jsonschema.exceptions import ValidationError

from .src import cfg, data, logs, modules, sge, util  # noqa: F401

_MP_MANAGER = mp.Manager()
_LOG_QUEUE = logs.get_log_queue(_MP_MANAGER)
_OUTPUT_QUEUE: mp.Queue = mp.Queue()
_PROCS: dict[UUID, modules.Runner] = {}
_STARTTIME = time.time()
_TIMESTAMP: str = time.strftime("%Y%m%d_%H%M%S", time.localtime(_STARTTIME))
CELLOPHANE_ROOT = Path(__file__).parent

click.rich_click.DEFAULT_STRING = "[{}]"


def _cleanup(logger):
    for proc in _PROCS.values():
        if proc.exitcode is None:
            logger.debug(f"Terminating {proc.label}")
            try:
                proc.terminate()
                proc.join()
            # Handle weird edge cases when terminating processes
            except Exception as exception:
                logger.debug(f"Failed to terminate {proc.label}: {exception}")


def _click_mapping(ctx, param, value):
    try:
        return cfg.parse_mapping(value)
    except Exception as exception:
        raise click.BadParameter(f"Invalid mapping: {exception}")


def _main(
    logger: logging.LoggerAdapter,
    config: cfg.Config,
    modules_path: Path,
    root: Path,
) -> None:
    """Run cellophane"""
    logger.setLevel(config.log_level)

    # Load modules
    hooks: list[type[modules.Hook]] = []
    runners: list[type[modules.Runner]] = []
    mixins: list[type[data.Mixin]] = []
    for base, obj in modules.load_modules(modules_path):
        if issubclass(obj, modules.Hook) and not obj == modules.Hook:
            logger.debug(f"Found hook {obj.__name__} ({base})")
            hooks.append(obj)
        elif issubclass(obj, data.Mixin) and not obj == data.Mixin:
            logger.debug(f"Found mixin {obj.__name__} ({base})")
            mixins.append(obj)
        elif issubclass(obj, modules.Runner) and not obj == modules.Runner:
            logger.debug(f"Found runner {obj.__name__} ({base})")
            runners.append(obj)

    # Resolve hook dependencies using topological sort
    try:
        hooks = modules.resolve_hook_dependencies(hooks)
    except Exception as exception:
        logger.error(f"Failed to resolve hook dependencies: {exception}")
        raise SystemExit(1)

    # Add mixins to data classes
    for mixin in [m for m in mixins]:
        logger.debug(f"Adding {mixin.__name__} mixin to samples")
        data.Samples.__bases__ = (*data.Samples.__bases__, mixin)
        if mixin.sample_mixin is not None:
            data.Sample.__bases__ = (*data.Sample.__bases__, mixin.sample_mixin)

    # Load samples from file, or create empty samples object
    if "samples_file" in config:
        samples = data.Samples.from_file(config.samples_file)
    else:
        samples = data.Samples()

    # Run pre-hooks
    for hook in [h() for h in hooks if h.when == "pre"]:
        result = hook(
            samples=deepcopy(samples),
            config=config,
            timestamp=_TIMESTAMP,
            log_queue=_LOG_QUEUE,
            log_level=config.log_level,
            root=root,
        )

        if issubclass(type(result), data.Samples):
            samples = result

    # Validate samples
    for invalid_sample in samples.validate():
        logger.warning(f"Removed invalid sample {invalid_sample.id}")

    # Start all loaded runners
    result_samples: data.Samples = data.Samples()
    sample_pids: dict[str, set[UUID]] = {s.id: set() for s in samples}
    try:
        for runner in runners if samples else []:
            logger.info(
                f"Starting runner {runner.__name__} for {len(samples)} samples"
            )

            for _samples in (
                samples.split(link_by=runner.link_by)
                if runner.individual_samples
                else [samples]
            ):
                proc = runner(
                    samples=_samples,
                    config=config,
                    timestamp=_TIMESTAMP,
                    output_queue=_OUTPUT_QUEUE,
                    log_queue=_LOG_QUEUE,
                    log_level=config.log_level,
                    root=root,
                )
                _PROCS[proc.id] = proc
                for sample in _samples:
                    sample_pids[sample.id] |= {proc.id}

        for proc in _PROCS.values():
            proc.start()

        # Wait for all runners to finish
        while not all(proc.done for proc in _PROCS.values()):
            result, pid = _OUTPUT_QUEUE.get()
            result_samples += result
            for sample in result:
                if sample.id in result_samples.complete.unique_ids:
                    logger.info(f"Sample {sample.id} completed by all runners")

            _PROCS[pid].join()
            _PROCS[pid].done = True

    except KeyboardInterrupt:
        logger.critical("Received SIGINT, telling runners to shut down...")
        _cleanup(logger)

    except Exception as e:
        logger.critical(f"Unhandled exception in runner: {e}")
        _cleanup(logger)

    finally:
        # Run post-hooks
        for hook in [h() for h in hooks if h.when == "post"]:
            match hook.condition:
                case "complete":
                    hook_samples = result_samples.complete
                case "failed":
                    hook_samples = result_samples.failed
                case "always":
                    hook_samples = result_samples
            if hook_samples:
                hook(
                    samples=hook_samples,
                    config=config,
                    timestamp=_TIMESTAMP,
                    log_queue=_LOG_QUEUE,
                    log_level=config.log_level,
                    root=root,
                )

        logger.info(
            f"Execution complete in {format_timespan(time.time() - _STARTTIME)}"
        )


def cellophane(
    label: str,
    wrapper_log: Optional[Path] = None,
    schema_path: Optional[Path] = None,
    modules_path: Optional[Path] = None,
) -> Callable:
    """Generate a cellophane CLI from a schema file"""
    root = Path(inspect.stack()[1].filename).parent
    _wrapper_log = wrapper_log or root / "pipeline.log"
    _schema_path = schema_path or root / "schema.yaml"
    _modules_path = modules_path or root / "modules"

    with (
        open(CELLOPHANE_ROOT / "schema.base.yaml", encoding="utf-8") as base_handle,
        open(_schema_path, "r", encoding="utf-8") as custom_handle,
    ):
        base = yaml.safe_load(base_handle)
        custom = yaml.safe_load(custom_handle)

    for module_schema_path in _modules_path.glob("*/schema.yaml"):
        with open(module_schema_path, "r", encoding="utf-8") as module_handle:
            module = yaml.safe_load(module_handle)
            custom = util.merge_mappings(custom, module)

    merged = util.merge_mappings(custom, base)
    schema = cfg.Schema(merged)

    @click.command()
    @logs.handle_logging(
        label=label,
        level=logging.INFO,
        path=_wrapper_log,
        queue=_LOG_QUEUE,
        propagate_exceptions=False,
    )
    def inner(config_path, logger, **kwargs) -> Any:
        _config = cfg.Config(config_path, schema, **kwargs)
        _config.analysis = label
        try:
            return _main(
                config=_config,
                logger=logger,
                modules_path=_modules_path,
                root=root,
            )
        except KeyboardInterrupt:
            logger.critical("Received SIGINT, telling active processes to shut down...")
            _cleanup(logger)
        except ValidationError as exception:
            _config = cfg.Config(config_path, schema, validate=False, **kwargs)
            for error in schema.iter_errors(_config):
                logger.critical(f"Invalid configuration: {error.message}")
            raise SystemExit(1) from exception

        except Exception as exception:
            logger.critical(
                f"Unhandled exception: {exception}",
                exc_info=_config.log_level == "DEBUG",
                stacklevel=2,
            )
            _cleanup(logger)

    for flag, _, default, description, secret, _type in schema.flags:
        inner = click.option(
            f"--{flag}",
            type=str if _type in (list, dict) else _type,
            callback=_click_mapping if _type == dict else None,
            is_flag=_type == bool,
            multiple=_type in (list, dict),
            default=default,
            help=description,
            show_default=not secret,
        )(inner)

    inner = click.option(
        "config_path",
        "--config",
        type=click.Path(exists=True),
        help="Path to config file",
        is_eager=True,
        callback=lambda ctx, _, value: cfg.set_defaults(ctx, value, schema),
    )(inner)

    return inner
