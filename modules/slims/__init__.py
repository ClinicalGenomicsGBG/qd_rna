"""Module for getting samples from SLIMS"""

from copy import deepcopy
from functools import cached_property, reduce
from json import loads
from logging import LoggerAdapter
from time import time
from typing import Optional, Any
from collections import UserList, UserDict
import re

from humanfriendly import parse_timespan
from slims.criteria import (
    Criterion,
    conjunction,
    disjunction,
    ends_with,
    starts_with,
    contains,
    is_one_of,
    equals,
    equals_ignore_case,
    greater_than,
    less_than,
    between_inclusive,
    is_not,
)
from slims.slims import Record, Slims

from cellophane import cfg, data, modules


def _parse_bool(string: str) -> list[str]:
    """
    Split string on "and"/"or" but not within parentheses

    >>> _parse_bool("a is x and (b is y or c is d) or g is h")
    ['a is x', 'and', 'b is y or c is d', 'or', 'g is h']
    """
    string = " ".join(string.split())
    parts: list[str] = []
    part = ""
    depth = 0
    while string:
        if string[0] == "(":
            if depth > 0:
                part += string[0]
            depth += 1
        elif string[0] == ")":
            if depth > 1:
                part += string[0]
            depth -= 1
        elif string[1:4] == "and" and depth == 0:
            parts.append(part)
            parts.append("and")
            string = string[4:]
            part = ""
        elif string[1:3] == "or" and depth == 0:
            parts.append(part)
            parts.append("or")
            string = string[3:]
            part = ""
        else:
            part += string[0]
            if len(string) == 1:
                parts.append(part)
        string = string[1:]

    while len(parts) == 1 and any(w in parts[0] for w in [" and ", " or ", "(", ")"]):
        parts = _parse_bool(parts[0])

    return parts


def _parse_criteria(criteria: str) -> Criterion:
    """Parse criteria"""

    _criteria = _parse_bool(criteria)
    if len(_criteria) == 1:
        _criteria = _criteria[0].split(" ")

    match _criteria:
        case [a, "and", b]:
            return conjunction().add(_parse_criteria(a)).add(_parse_criteria(b))
        case [a, "or", b]:
            return disjunction().add(_parse_criteria(a)).add(_parse_criteria(b))
        case [field, *_] if not field.startswith("cntn_"):
            raise ValueError(f"Invalid field: {field}")

        case [field, "equals", value]:
            return equals(field, value)
        case [field, "not_equals", value]:
            return is_not(equals(field, value))
        case [field, "is_one_of", *values]:
            return is_one_of(field, values)

        case [field, "equals_ignore_case", value]:
            return equals_ignore_case(field, value)
        case [field, "not_equals_ignore_case ", value]:
            return is_not(equals_ignore_case(field, value))

        case [field, "contains", value]:
            return contains(field, value)
        case [field, "not_contains ", value]:
            return is_not(contains(field, value))

        case [field, "starts_with", value]:
            return starts_with(field, value)
        case [field, "not_starts_with", value]:
            return is_not(starts_with(field, value))

        case [field, "ends_with", value]:
            return ends_with(field, value)
        case [field, "not_ends_with", value]:
            return is_not(ends_with(field, value))

        case [field, "between_inclusive", *values]:
            return between_inclusive(field, *values)
        case [field, "not_between_inclusive", *values]:
            return is_not(between_inclusive(field, *values))

        case [field, "greater_than", value]:
            return greater_than(field, value)
        case [field, "less_than", value]:
            return less_than(field, value)
        case _:
            raise ValueError(f"Invalid criteria: {criteria}")


def _get_field(record: Record, field: str, default=None) -> Any:
    """Get a field from SLIMS record"""
    try:
        if field.startswith("json:"):
            _field, *_key = re.split(r"\.|(\[[0-9]*\])", field[5:])
            _key = [int(k.strip("[]")) if k.startswith("[") else k for k in _key if k]
            _json = loads(record.__dict__[_field].value)
            return reduce(lambda x, y: x[y], _key, _json)
        else:
            return record.__getattribute__(field).value
    except AttributeError:
        return default
    except KeyError:
        return default


def get_records(
    *args: Criterion,
    connection: Slims,
    slims_id: Optional[str | list[str]] = None,
    max_age: Optional[int | str] = None,
    content_type: Optional[int | list[int]] = None,
    derived_from: Optional[str | list[str]] = None,
    **kwargs: str | int | list[str | int],
) -> list[Record]:
    """Get records from SLIMS"""

    criteria = conjunction()

    match slims_id:
        case str():
            criteria = criteria.add(equals("cntn_id", slims_id))
        case [*ids]:
            criteria = criteria.add(is_one_of("cntn_id", ids))
        case _ if slims_id is not None:
            raise TypeError(f"Invalid type for id: {type(slims_id)}")

    match max_age:
        case int() | str():
            min_mtime = int(time() - parse_timespan(str(max_age))) * 1e3
            criteria = criteria.add(greater_than("cntn_modifiedOn", min_mtime))
        case _ if max_age is not None:
            raise TypeError(f"Expected int or str, got {type(max_age)}")

    match content_type:
        case None:
            pass
        case int():
            criteria = criteria.add(equals("cntn_fk_contentType", content_type))
        case [*_, int()]:
            criteria = criteria.add(is_one_of("cntn_fk_contentType", content_type))
        case _:
            raise TypeError(f"Expected int(s), got {type(content_type)}")

    match derived_from:
        case record if isinstance(record, Record):
            original = {record.pk(): record}
            criteria = criteria.add(is_one_of("cntn_fk_originalContent", [*original]))
        case [*records] if all(isinstance(r, Record) for r in records):
            original = {r.pk(): r for r in records}
            criteria = criteria.add(is_one_of("cntn_fk_originalContent", [*original]))
        case _:
            raise TypeError(f"Expected Record(s), got {derived_from}")

    for key, value in kwargs.items():
        criteria = criteria.add(
            is_one_of(key, [value] if isinstance(value, int | str) else value)
        )

    for arg in args:
        criteria = criteria.add(arg)

    return connection.fetch("Content", criteria)


class SlimsSample:
    """A sample container with SLIMS integration"""

    record: Optional[Record]
    bioinformatics: Optional[Record]

    @classmethod
    def from_record(cls, record: Record, config: cfg.Config, **kwargs):
        """Create a sample from a SLIMS fastq record"""
        return cls(
            id=record.cntn_id.value,
            bioinformatics=None,
            record=record,
            **{
                key: _get_field(record, field)
                for key, field in config.slims.field.items()
            },
            **kwargs,
        )

    def add_bioinformatics(
        self,
        config: cfg.Config,
    ):
        """Add a bioinformatics record to the sample"""

        if self.bioinformatics is None and self._connection is not None:
            fields = {
                "cntn_id": self.record.cntn_id.value,
                "cntn_fk_contentType": config.slims.bioinfo.content_type,
                "cntn_status": 10,  # Pending
                "cntn_fk_location": 83,  # FIXME: Should location be configuarable?
                "cntn_fk_originalContent": self.record.pk(),
                "cntn_fk_user": "",  # FIXME: Should user be configuarable?
                config.slims.bioinfo.state_field: "novel",
            }

            self.bioinformatics = self._connection.add("Content", fields)

    def set_bioinformatics_state(self, state, config: cfg.Config):
        """Set the bioinformatics state"""

        match state:
            case "running" | "complete" | "error":
                if self.bioinformatics is not None:
                    self.bioinformatics = self.bioinformatics.update(
                        {config.slims.bioinfo.state_field: state}
                    )
            case _:
                raise ValueError(f"Invalid state: {state}")

    @cached_property
    def _connection(self) -> Optional[Slims]:
        """Get a connection to SLIMS from the record"""

        return (
            Slims(
                "cellophane",
                url=self.record.slims_api.raw_url,
                username=self.record.slims_api.username,
                password=self.record.slims_api.password,
            )
            if self.record is not None
            else None
        )

    def __reduce__(self) -> str | tuple:
        """Remove open connection before pickle"""

        if hasattr(self, "_connection"):
            delattr(self, "_connection")
        return super().__reduce__()


class SlimsSamples(data.Mixin, sample_mixin=SlimsSample):
    """A list of sample containers with SLIMS integration"""

    @classmethod
    def from_records(
        cls,
        records: list[Record],
        config: cfg.Config,
    ) -> "SlimsSamples":
        """Get samples from SLIMS records"""

        return cls(
            [
                cls.sample_class.from_record(record=record, config=config)
                for record in records
            ]
        )

    def add_bioinformatics(self, config: cfg.Config) -> None:
        """Add bioinformatics content to SLIMS samples"""
        for sample in self:
            sample.add_bioinformatics(config)

    def set_bioinformatics_state(self, state: str, config: cfg.Config) -> None:
        """Update bioinformatics state in SLIMS"""
        match state:
            case "running" | "complete" | "error":
                for sample in self:
                    sample.set_bioinformatics_state(state, config)
            case _:
                raise ValueError(f"Invalid state: {state}")


@modules.pre_hook(label="SLIMS Fetch", priority=0)
def slims_samples(
    samples: data.Samples,
    config: cfg.Config,
    logger: LoggerAdapter,
    **_,
) -> Optional[SlimsSamples]:
    """Load novel samples from SLIMS."""
    if "slims" in config:
        slims_connection = Slims(
            name=__package__,
            url=config.slims.url,
            username=config.slims.username,
            password=config.slims.password,
        )
        
        if config.slims.derived_from:
            logger.debug("Fetching derived from records")
            parent_records = get_records(
                _parse_criteria(config.slims.derived_from.criteria),
                content_type=config.slims.derived_from.content_type,
                connection=slims_connection,
            )
            samples_kwargs = {"derived_from": parent_records}


        if samples:
            logger.debug("Augmenting existing samples")
            samples_kwargs = {"slms_id": [s.id for s in samples]}

        elif config.slims.get("id", None):
            logger.debug("Fetching samples by ID")
            samples_kwargs = {"slms_id": config.slims.id}

        else:
            logger.debug(f"Fetching samples from the last {config.slims.novel_max_age}")
            samples_kwargs = {"max_age": config.slims.novel_max_age}


        logger.debug(f"Fetching samples with {samples_kwargs}")
        records = get_records(
            _parse_criteria(config.slims.criteria),
            content_type=config.slims.content_type,
            connection=slims_connection,
            **samples_kwargs,
        )


        if config.slims.bioinfo.check:
            bioinfo = get_records(
                _parse_criteria(config.slims.bioinfo.check_criteria),
                connection=slims_connection,
                derived_from=records,
                content_type=config.slims.bioinfo.content_type,
            )

            records = [
                record
                for record in records
                if record.pk() not in [
                    b.cntn_fk_originalContent.value
                    for b in bioinfo
                ]
            ]
            
            logger.debug(f"skipping {len(bioinfo)} samples with complete bioinformatics")

        slims_samples = samples.from_records(records, config)

        if samples:
            for idx, sample in enumerate(samples):
                match = [m for m in slims_samples if m.id == sample.id]
                common_keys = set([k for s in match for k in s]) & set(sample.keys())
                common_keys -= set(["files", "backup"])
                for key in common_keys:
                    _match = []
                    for match_sample in match:
                        if (
                            (m_value := match_sample[key])
                            and (s_value := sample[key])
                            and s_value == m_value
                        ):
                            _match.append(match_sample)
                    match = _match

                if len(match) > 1:
                    logger.warning(f"Multiple SLIMS samples found for {sample.id}")
                elif len(match) == 0:
                    logger.warning(f"SLIMS sample not found for {sample.id}")
                else:
                    if sample.files == None:
                        sample.pop("files")
                    _data = {**match[0]} | {**sample}

                    samples[idx] = sample.__class__(
                        id=_data.pop("id"),
                        **deepcopy(_data),
                    )
            return samples
        else:
            return slims_samples

    else:
        logger.warning("No SLIMS connection configured")
        return None


@modules.pre_hook(label="SLIMS Add Bioinfo")
def slims_bioinformatics(
    samples: data.Samples,
    config: cfg.Config,
    logger: LoggerAdapter,
    **_,
) -> None:
    """Load novel samples from SLIMS."""
    if config.slims.dry_run:
        logger.debug("Dry run - Not adding bioinformatics")
    elif config.slims.bioinfo.create:
        logger.info("Adding bioinformatics content")
        samples.add_bioinformatics(config)
        samples.set_bioinformatics_state("running", config)
    return samples


@modules.post_hook(label="SLIMS Update Bioinfo")
def slims_update(
    config: cfg.Config,
    samples: data.Samples,
    logger: LoggerAdapter,
    **_,
) -> None:
    """Update SLIMS samples with bioinformatics content."""

    if config.slims.dry_run:
        logger.info("Dry run - Not updating SLIMS")
    elif isinstance(samples, SlimsSamples):
        logger.info("Updating bioinformatics")
        unique = {
            pk: [s for s in samples if s.record.pk() == pk]
            for pk in set(s.record.pk() for s in samples)
        }
        for _samples in unique.values():
            if all(s.complete for s in _samples):
                _samples[0].set_bioinformatics_state("complete", config)
            else:
                _samples[0].set_bioinformatics_state("error", config)
    else:
        logger.info("No SLIMS samples to update")
