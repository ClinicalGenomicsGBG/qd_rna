"""Cellophane: A library for writing modular wrappers"""
import inspect
import logging
import multiprocessing as mp
import sys
import time
from copy import deepcopy
from graphlib import CycleError, TopologicalSorter
from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path
from typing import Any, Callable, Iterator, Optional
from uuid import UUID

import rich_click as click
import yaml
from humanfriendly import format_timespan
from jsonschema.exceptions import ValidationError

from .src import cfg, data, logs, modules, sge, util

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


def _load_modules(
    path: Path,
) -> Iterator[tuple[str, type[modules.Hook] | type[modules.Runner] | type[data.Mixin]]]:
    for file in [*path.glob("*.py"), *path.glob("*/__init__.py")]:
        base = file.stem if file.stem != "__init__" else file.parent.name
        name = f"_cellophane_module_{base}"
        spec = spec_from_file_location(name, file)
        original_handlers = logging.root.handlers.copy()
        if spec is not None:
            module = module_from_spec(spec)
            if spec.loader is not None:
                try:
                    sys.modules[name] = module
                    spec.loader.exec_module(module)
                except ImportError:
                    pass
                else:
                    # Reset logging handlers to avoid duplicate messages
                    for handler in logging.root.handlers:
                        if handler not in original_handlers:
                            handler.close()
                            logging.root.removeHandler(handler)

                    for obj in [getattr(module, a) for a in dir(module)]:
                        if isinstance(obj, type) and (
                            issubclass(obj, modules.Hook)
                            or issubclass(obj, data.Mixin)
                            or issubclass(obj, modules.Runner)
                        ) and inspect.getmodule(obj) == module:
                            yield base, obj


def _resolve_hook_dependencies(
    hooks: list[type[modules.Hook]],
) -> list[type[modules.Hook]]:
    deps = {
        name: {
            *[d for h in hooks if h.__name__ == name for d in h.after],
            *[h.__name__ for h in hooks if name in h.before],
        }
        for name in {
            *[n for h in hooks for n in h.before + h.after],
            *[h.__name__ for h in hooks],
        }
    }

    order = [*TopologicalSorter(deps).static_order()]
    return [*sorted(hooks, key=lambda h: order.index(h.__name__))]


def _main(
    logger: logging.LoggerAdapter,
    config: cfg.Config,
    modules_path: Path,
    root: Path,
) -> None:
    """Run cellophane"""
    logger.setLevel(config.log_level)

    hooks: list[type[modules.Hook]] = []
    runners: list[type[modules.Runner]] = []
    mixins: list[type[data.Mixin]] = []
    for base, obj in _load_modules(modules_path):
        if issubclass(obj, modules.Hook) and not obj == modules.Hook:
            logger.debug(f"Found hook {obj.__name__} ({base})")
            hooks.append(obj)
        elif issubclass(obj, data.Mixin) and not obj == data.Mixin:
            logger.debug(f"Found mixin {obj.__name__} ({base})")
            mixins.append(obj)
        elif issubclass(obj, modules.Runner) and not obj == modules.Runner:
            logger.debug(f"Found runner {obj.__name__} ({base})")
            runners.append(obj)

    hooks = _resolve_hook_dependencies(hooks)

    for mixin in [m for m in mixins]:
        logger.debug(f"Adding {mixin.__name__} mixin to samples")
        data.Samples.__bases__ = (*data.Samples.__bases__, mixin)
        if mixin.sample_mixin is not None:
            data.Sample.__bases__ = (*data.Sample.__bases__, mixin.sample_mixin)
    
    if "samples_file" in config:
        samples = data.Samples.from_file(config.samples_file)
    else:
        samples = data.Samples()
    
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

    for invalid_sample in samples.validate():
        logger.warning(f"Removed invalid sample {invalid_sample.id}")

    
    result_samples = data.Samples()
    sample_pids: dict(str, set[UUID]) = {s.id: set() for s in samples}

    failed_samples = data.Samples()
    partial_samples = data.Samples()
    complete_samples = data.Samples()
    try:
        if samples:
            for runner in runners:
                logger.info(f"Starting runner {runner.__name__} for {len(samples)} samples")

                for _samples in samples.split(link_by=runner.link_by) if runner.individual_samples else [samples]:
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

        while not all(proc.done for proc in _PROCS.values()):
            result, pid = _OUTPUT_QUEUE.get()
            result_samples += result
            for sample in result:
                if sample.done:
                    sample_pids[sample.id] -= {pid}
                    if not sample_pids[sample.id]:
                        logger.info(f"Sample {sample.id} completed by all runners")
                        complete_samples += [s for s in result_samples if s.id == sample.id]
                        partial_samples = data.Samples(s for s in result_samples if s.id == sample.id)
                        
                    else:
                        partial_samples += [sample]
                else:
                    failed_samples += [sample]

            _PROCS[pid].join()
            _PROCS[pid].done = True
        
    except KeyboardInterrupt:
        logger.critical("Received SIGINT, telling runners to shut down...")
        _cleanup(logger)
    
    except Exception as e:
        logger.critical(f"Unhandled exception in runner: {e}")
        _cleanup(logger)

    finally:
        failed_samples += data.Samples(s for s in samples if s.id not in [r.id for r in result_samples])
        for hook in [h() for h in hooks if h.when == "post"]:

            hook(
                samples=data.Samples(
                    [
                        *complete_samples,
                        *(partial_samples if hook.condition != "complete" else []),
                        *(failed_samples if hook.condition == "always" else []),
                    ]
                ),
                config=config,
                timestamp=_TIMESTAMP,
                log_queue=_LOG_QUEUE,
                log_level=config.log_level,
                root=root,
            )

        logger.info(f"Execution complete in {format_timespan(time.time() - _STARTTIME)}")


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
        try:
            return _main(
                config=_config,
                logger=logger,
                modules_path=_modules_path,
                root=root,
            )
        except CycleError as exception:
            logger.error(f"Circular dependency in hooks: {exception}")
            raise SystemExit(1)

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
