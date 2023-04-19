"""Configuration file handling"""

from copy import deepcopy
from functools import reduce
from pathlib import Path
from typing import Sequence, Optional, Type, Iterator, Mapping, Any

import rich_click as click
from jsonschema import Draft7Validator, Validator, validators, ValidationError
from yaml import safe_load

from . import data, util


def _set_options(cls, validate: bool) -> Type[Validator]:
    """Set default values when validating"""
    validate_properties = cls.VALIDATORS["properties"]

    def _is_object_or_container(_, instance):
        return cls.TYPE_CHECKER.is_type(instance, "object") or isinstance(
            instance, data.Container
        )

    def _is_array(_, instance):
        return cls.TYPE_CHECKER.is_type(instance, "array") or isinstance(
            instance, Sequence
        )

    def _is_path(_, instance):
        return isinstance(instance, Optional[Path | click.Path])

    def _is_mapping(_, instance):
        return isinstance(instance, Sequence) and all(
            isinstance(i, Mapping) for i in instance
        )

    def _set(cls, properties, instance, schema):
        for prop, subschema in properties.items():
            if "default" in subschema:
                instance.setdefault(prop, subschema["default"])
            if subschema.get("type", None) == "path" and prop in instance:
                instance[prop] = Path(instance[prop])

        if validate:
            for error in validate_properties(cls, properties, instance, schema):
                yield error

    type_checker = cls.TYPE_CHECKER.redefine_many(
        {
            "object": _is_object_or_container,
            "array": _is_array,
            "path": _is_path,
            "mapping": _is_mapping,
        }
    )
    return validators.extend(cls, {"properties": _set}, type_checker=type_checker)


def _get_options(cls) -> Type[Validator]:
    """Parse options from schema and ignore validation errors"""
    validate_properties = cls.VALIDATORS["properties"]

    def func(cls, properties, instance, schema):
        for prop, subschema in properties.items():
            if "properties" in subschema:
                instance[prop] = {}
            else:
                instance[prop] = []
                match subschema:
                    case {"default": default}:
                        instance[prop].append(default)
                    case _:
                        instance[prop].append(None)

                match subschema:
                    case {"description": description}:
                        instance[prop].append(description)
                    case _:
                        instance[prop].append("")

                match subschema:
                    case {"secret": True}:
                        instance[prop].append(True)
                    case _:
                        instance[prop].append(False)

                match subschema:
                    case {"enum": enum}:
                        instance[prop].append(click.Choice(enum, case_sensitive=False))
                    case {"type": "boolean"}:
                        instance[prop].append(bool)
                    case {"type": "path"}:
                        instance[prop].append(click.Path())
                    case {"type": "string"}:
                        instance[prop].append(str)
                    case {"type": "integer"}:
                        instance[prop].append(int)
                    case {"type": "number"}:
                        instance[prop].append(float)
                    case {"type": "array"}:
                        instance[prop].append(list)
                    case {"type": "mapping"}:
                        instance[prop].append(dict)
                    case _:
                        instance[prop].append(None)

        for error in validate_properties(cls, properties, instance, schema):
            yield error

    return validators.extend(
        cls, {k: None for k in cls.VALIDATORS} | {"properties": func}
    )


def parse_mapping(string_mapping: dict | Sequence[str] | str) -> list[dict[str, Any]]:
    match string_mapping:
        case dict(mapping):
            return [mapping]
        case [*strings]:
            return [m for s in strings for m in parse_mapping(s)]
        case str(string):
            _mapping: dict[str, Any] = {}
            for kv in string.split():
                for k, v in [kv.split("=")]:
                    identifier = k.strip("{}")
                    if not identifier.isidentifier():
                        raise ValueError(f"{identifier} is not a valid identifier")
                    else:
                        _mapping[identifier] = v
            return [_mapping]
        case _:
            raise ValueError("format must be 'key=value ...'")


class Schema(data.Container):
    """Schema for validating configuration files"""

    @classmethod
    def from_file(cls, path: Path | Sequence[Path]):
        """Load schema from file"""
        path = [path] if isinstance(path, Path) else path
        schema: dict = {}
        for file in path:
            with open(file, "r", encoding="utf-8") as handle:
                schema = util.merge_mappings(
                    deepcopy(schema),
                    safe_load(handle) or {},
                )
        cls(schema)

    @property
    def properties(self) -> dict:
        """Get properties from schema"""
        _properties: dict = {}
        _get_options(Draft7Validator)({**self.data}).validate(_properties)
        return _properties

    @property
    def key_map(self) -> list[list[str]]:
        """Get key map from schema"""
        return util.map_nested_keys(self.properties)

    @property
    def flags(self):
        """Get flags from schema"""
        for key in self.key_map:  # pylint: disable=not-an-iterable
            flag = "_".join(key)
            default, description, secret, _type = reduce(
                lambda x, y: x[y], key, self.properties
            )
            yield flag, key, default, description, secret, _type

    def set_defaults(
        self,
        config: data.Container,
        validate: bool = True,
    ) -> data.Container:
        """Validate configuration"""
        _config = deepcopy(config)
        _validator = _set_options(Draft7Validator, validate=validate)({**self.data})
        _validator.validate(_config)
        return _config

    def iter_errors(self, config: data.Container) -> Iterator[ValidationError]:
        """Iterate over validation errors"""
        _config = deepcopy(config)
        _validator = _set_options(Draft7Validator, validate=True)({**self.data})
        return _validator.iter_errors(_config)


class Config(data.Container):
    """Configuration file"""

    def __init__(
        self,
        path: Optional[Path | click.Path],
        schema: Schema,
        validate: bool = True,
        **kwargs,
    ):

        if path is not None:
            with open(str(path), "r", encoding="utf-8") as handle:
                _data = data.Container(safe_load(handle))
        else:
            _data = data.Container({})

        for flag, key, *_ in schema.flags:
            if flag not in _data and kwargs[flag] is not None:
                _data[key] = kwargs[flag]

        _data = schema.set_defaults(_data, validate=validate)
        super().__init__(_data)


def set_defaults(
    ctx: click.Context,
    config_path: Optional[Path | click.Path],
    schema: Schema,
):
    """Set click defaults from schema and configuration file"""
    if config_path is not None:
        with open(str(config_path), "r", encoding="utf-8") as handle:
            config = data.Container(safe_load(handle))
    else:
        config = data.Container({})

    ctx.default_map = {
        flag: config[key] if key in config else default
        for flag, key, default, *_ in schema.flags
    }
