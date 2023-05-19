"""Module for getting samples from SLIMS"""

import re
from copy import deepcopy
from datetime import datetime, timedelta
from functools import cached_property, reduce
from json import loads
from logging import LoggerAdapter
from typing import Any, Literal

from attrs import define, field
from cellophane import cfg, data, modules
from humanfriendly import parse_timespan

from slims.criteria import (
    Criterion,
    between_inclusive,
    conjunction,
    contains,
    disjunction,
    ends_with,
    equals,
    equals_ignore_case,
    greater_than,
    is_not,
    is_one_of,
    less_than,
    starts_with,
)
from slims.slims import Record, Slims


def _split_criteria(criteria: str) -> list[str]:
    """
    Split string on "and"/"or" but not within parentheses

    >>> _parse_bool("a is x and (b is y or c is d) or g is h")
    ['a is x', 'and', 'b is y or c is d', 'or', 'g is h']
    """

    _criteria = " ".join(criteria.split())
    parts: list[str] = []
    part = ""
    depth = 0
    while _criteria:
        if _criteria[0] == "(":
            if depth > 0:
                part += "("
            depth += 1
        elif _criteria[0] == ")":
            if depth > 1:
                part += ")"
            depth -= 1
            parts += (part,)
        elif _criteria[1:4] == "and" and depth == 0:
            parts += part, "and"
            _criteria, part = _criteria[4:], ""
        elif _criteria[1:3] == "or" and depth == 0:
            parts += part, "or"
            _criteria, part = _criteria[3:], ""
        else:
            part += _criteria[0]
            if len(_criteria) == 1:
                parts.append(part)
        _criteria = _criteria[1:]

    if depth != 0:
        raise ValueError(f"Unmatched parentheses: {criteria}")

    while len(parts) == 1 and any(w in parts[0] for w in [" and ", " or ", "(", ")"]):
        parts = _split_criteria(parts[0])

    return parts


def _parse_criteria(  # type: ignore[return]
    criteria: str | list[str],
    parent_records: list[Record] | None = None,
) -> list[Criterion]:
    """Parse criteria"""

    match criteria:
        case str(criteria) if criteria.startswith("->") and not parent_records:
            raise ValueError("Cannot use leading '->' without parent record(s)")
        case str(criteria) if criteria.startswith("->"):
            _parsed = _parse_criteria(criteria[2:])
            if parent_records is None:
                raise ValueError("Cannot use leading '->' without parent record(s)")
            else:
                parent_pks = [p.pk() for p in parent_records]
                return [
                    conjunction()
                    .add(is_one_of("cntn_fk_originalContent", parent_pks))
                    .add(_parsed[0]),
                    *_parsed[1:],
                ]
        case str(criteria) if parent_records:
            _parsed = _parse_criteria(criteria)
            parent_pks = [p.pk() for p in parent_records]
            return [
                conjunction().add(is_one_of("cntn_pk", parent_pks)).add(_parsed[0]),
                *_parsed[1:],
            ]
        case str(criteria) if "->" in criteria:
            return [_parse_criteria(c)[0] for c in criteria.split("->")]
        case str(criteria):
            return _parse_criteria(_split_criteria(criteria))

        case [criterion]:
            criteria = _split_criteria(criterion)
            if len(criteria) == 1:
                return _parse_criteria(criteria[0].split(" "))
            else:
                return _parse_criteria(criteria)

        case [a, "and", *b]:
            return [conjunction().add(_parse_criteria(a)[0]).add(_parse_criteria(b)[0])]
        case [a, "or", *b]:
            return [disjunction().add(_parse_criteria(a)[0]).add(_parse_criteria(b)[0])]
        case [field, *_] if not field.startswith("cntn_"):
            raise ValueError(f"Invalid field: {field}")

        case [field, "equals", value]:
            return [equals(field, value)]
        case [field, "not_equals", value]:
            return [is_not(equals(field, value))]
        case [field, "one_of", *values]:
            return [is_one_of(field, values)]
        case [field, "not_one_of", *values]:
            return [is_not(is_one_of(field, values))]

        case [field, "equals_ignore_case", value]:
            return [equals_ignore_case(field, value)]
        case [field, "not_equals_ignore_case ", value]:
            return [is_not(equals_ignore_case(field, value))]

        case [field, "contains", value]:
            return [contains(field, value)]
        case [field, "not_contains ", value]:
            return [is_not(contains(field, value))]

        case [field, "starts_with", value]:
            return [starts_with(field, value)]
        case [field, "not_starts_with", value]:
            return [is_not(starts_with(field, value))]

        case [field, "ends_with", value]:
            return [ends_with(field, value)]
        case [field, "not_ends_with", value]:
            return [is_not(ends_with(field, value))]

        case [field, "between", *values]:
            return [between_inclusive(field, *values)]
        case [field, "not_between", *values]:
            return [is_not(between_inclusive(field, *values))]

        case [field, "greater_than", value]:
            return [greater_than(field, value)]
        case [field, "less_than", value]:
            return [less_than(field, value)]
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
    *args,
    connection: Slims,
    string_criteria: str | None = None,
    slims_id: str | list[str] | None = None,
    content_type: int | list[int] | None = None,
    max_age: int | str | None = None,
    derived_from: Record | list[Record] | None = None,
    unrestrict_parents: bool = False,
    **kwargs: str | int | list[str | int],
) -> list[Record]:
    """Get records from SLIMS"""

    if isinstance(derived_from, Record):
        derived_from = [derived_from]

    if string_criteria:
        _parsed = _parse_criteria(string_criteria, parent_records=derived_from)
        parent_records: list[Record] | None = None
        for criterion in _parsed[:-1]:
            parent_records = get_records(
                criterion,
                connection=connection,
                derived_from=parent_records,
                slims_id=slims_id if not unrestrict_parents else None,
                max_age=max_age if not unrestrict_parents else None,
            )

        return get_records(
            _parsed[-1],
            connection=connection,
            derived_from=parent_records,
            slims_id=slims_id,
            max_age=max_age,
        )

    else:
        criteria = conjunction()
        match slims_id:
            case str():
                criteria = criteria.add(equals("cntn_id", slims_id))
            case [*ids]:
                criteria = criteria.add(is_one_of("cntn_id", ids))
            case _ if slims_id is not None:
                raise TypeError(f"Invalid type for id: {type(slims_id)}")

        match content_type:
            case int():
                criteria = criteria.add(equals("cntn_fk_contentType", content_type))
            case [*types]:
                criteria = criteria.add(is_one_of("cntn_fk_contentType", types))
            case _ if content_type is not None:
                raise TypeError(f"Invalid type for content_type: {type(content_type)}")

        match max_age:
            case int() | str():
                now = datetime.now()
                max_date = now - timedelta(seconds=parse_timespan(str(max_age)))
                criteria = criteria.add(
                    between_inclusive("cntn_modifiedOn", max_date, now)
                )
            case _ if max_age is not None:
                raise TypeError(f"Expected int or str, got {type(max_age)}")

        match derived_from:
            case None:
                pass
            case record if isinstance(record, Record):
                original = {record.pk(): record}
                criteria = criteria.add(
                    is_one_of("cntn_fk_originalContent", [*original])
                )
            case [*records] if all(isinstance(r, Record) for r in records):
                original = {r.pk(): r for r in records}
                criteria = criteria.add(
                    is_one_of("cntn_fk_originalContent", [*original])
                )
            case _:
                raise TypeError(f"Expected Record(s), got {derived_from}")

        for key, value in kwargs.items():
            criteria = criteria.add(
                is_one_of(key, [value] if isinstance(value, int | str) else value)
            )

        for arg in args:
            criteria = criteria.add(arg)

        return connection.fetch("Content", criteria)


@define(slots=False, init=False)
class SlimsSample(data.Sample):
    """A sample container with SLIMS integration"""

    derived: list[tuple[Record, dict]] | None = field(default=None)
    record: Record | None = field(default=None)
    state: Literal["novel", "running", "complete", "error"] = field(default="novel")

    @classmethod
    def from_record(cls, record: Record, config: cfg.Config, **kwargs):
        """Create a sample from a SLIMS fastq record"""

        _sample = cls(
            id=record.cntn_id.value,
            state="novel",
            **{
                key: _get_field(record, field)
                for key, field in config.slims.map[0].items()
            },
            **kwargs,
        )  # type: ignore[call-arg]
        _sample.record = record
        return _sample

    def update_derived(
        self,
        config: cfg.Config,
    ):
        """Update/add derived records for the sample"""
        if not self.derived:
            self.derived = [(None, key_map) for key_map in config.slims.derive]
        if self.record:
            for idx, (record, key_map) in enumerate(self.derived):
                fields = {key: value.format(**self) for key, value in key_map.items()}
                fields |= {
                    "cntn_id": self.record.cntn_id.value,
                    "cntn_fk_originalContent": self.record.pk(),
                    "cntn_fk_user": config.slims.username,
                }
                if record:
                    self.derived[idx] = (record.update(fields), key_map)
                elif self._connection:
                    self.derived[idx] = (
                        self._connection.add("Content", fields),
                        key_map,
                    )

    @derived.validator
    def validate_derived(
        self,
        attribute: str,
        value: list[tuple[Record | None, dict]] | None,
    ):
        if not (value is None or isinstance(value, list)):
            raise ValueError(f"Expected 'None|list', got {value}")
        elif value is not None and not all(
            isinstance(v, tuple)
            and len(v) == 2
            and (isinstance(v[0], Record) or v[0] is None)
            and isinstance(v[1], dict)
            for v in value
        ):
            raise ValueError(
                f"Expected 'list[tuple[Record|None, dict]' for {attribute}, got {value}"
            )

    @record.validator
    def validate_record(self, attribute: str, value: Record | None):
        if not (value is None or isinstance(value, Record)):
            raise ValueError(
                f"Expected 'NoneType' or 'Record' for {attribute}, got {value}"
            )

    @state.validator
    def validate_state(
        self, attribute: str, value: Literal["novel", "running", "complete", "error"]
    ):
        """Set the state of the sample"""
        if value not in ["novel", "running", "complete", "error"]:
            raise ValueError(f"Invalid value for {attribute}: {value}")

    @cached_property
    def _connection(self) -> Slims | None:
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


class SlimsSamples(data.Samples):
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

    def update_derived(
        self,
        config: cfg.Config,
    ) -> None:
        """Update derived records in SLIMS"""
        for sample in self:
            sample.update_derived(config)

    def set_state(self, value: Literal["novel", "running", "complete", "error"]):
        """Set the state of the samples"""
        if value not in ["novel", "running", "complete", "error"]:
            raise ValueError(f"Invalid state: {value}")
        else:
            for sample in self:
                sample.state = value


@modules.pre_hook(label="SLIMS Fetch", before=["hcp_fetch", "slims_bioinformatics"])
def slims_fetch(
    samples: data.Samples,
    config: cfg.Config,
    logger: LoggerAdapter,
    **_,
) -> SlimsSamples | None:
    """Load novel samples from SLIMS."""
    if "slims" in config:
        slims_connection = Slims(
            name=__package__,
            url=config.slims.url,
            username=config.slims.username,
            password=config.slims.password,
        )

        slims_ids: list[str] | None = None
        max_age: str | None = None
        if samples:
            logger.info("Augmenting existing samples with info from SLIMS")
            slims_ids = [s.id for s in samples]

        elif config.slims.id:
            logger.info("Fetching samples from SLIMS by ID")
            slims_ids = config.slims.id

        else:
            logger.info(f"Fetching samples from the last {config.slims.novel_max_age}")
            max_age = config.slims.novel_max_age

        records = get_records(
            string_criteria=config.slims.find_criteria,
            connection=slims_connection,
            slims_id=slims_ids,
            max_age=max_age,
            unrestrict_parents=config.slims.unrestrict_parents,
        )

        if samples and records:
            for idx, sample in enumerate(samples):
                match = [
                    m
                    for m in samples.from_records(records, config)
                    if m.id == sample.id
                ]
                common_keys = set([k for s in match for k in s]) & set(sample.keys())
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
                    if sample.files is None:
                        sample.pop("files")
                    _data = {**match[0]} | {**sample}

                    samples[idx] = sample.__class__(
                        id=_data.pop("id"),
                        **deepcopy(_data),
                    )

            return samples

        elif records:
            if config.slims.check:
                logger.info("Checking SLIMS for completed samples")
                check = get_records(
                    string_criteria=config.slims.check_criteria,
                    connection=slims_connection,
                    derived_from=records,
                )

                original_ids = [r.cntn_id.value for r in records]
                records = [
                    record
                    for record in records
                    if record.pk()
                    not in [b.cntn_fk_originalContent.value for b in check]
                ]

                for sid in set(original_ids) - set([r.cntn_id.value for r in records]):
                    logger.info(f"Found completed bioinformatics for {sid}")

            return samples.from_records(records, config)

        else:
            logger.warning("No SLIMS samples found")
            return None

    else:
        logger.warning("No SLIMS connection configured")
        return None


@modules.pre_hook(label="SLIMS Derive", after=["slims_fetch"])
def slims_derive(
    samples: SlimsSamples,
    config: cfg.Config,
    logger: LoggerAdapter,
    **_,
) -> SlimsSamples:
    """Add derived content to SLIMS samples"""
    if config.slims.dry_run:
        logger.debug("Dry run - Not adding derived records")
    elif config.slims.derive and samples:
        logger.info("Creating derived records")
        samples.update_derived(config)
    return samples


@modules.pre_hook(label="SLIMS Mark Running", after="all")
def slims_running(
    samples: SlimsSamples,
    config: cfg.Config,
    logger: LoggerAdapter,
    **_,
) -> SlimsSamples:
    """Add derived content to SLIMS samples"""
    samples.set_state("running")
    if not config.slims.dry_run:
        logger.info("Setting SLIMS samples to running")
        samples.update_derived(config)
    return samples


@modules.post_hook(label="SLIMS Update Derived")
def slims_update(
    config: cfg.Config,
    samples: SlimsSamples,
    logger: LoggerAdapter,
    **_,
) -> None:
    """Update SLIMS samples and derived records."""
    if config.slims.dry_run:
        logger.info("Dry run - Not updating SLIMS")
        return

    complete = samples.__class__.from_records(
        records=[*set(s.record for s in samples.complete)], config=config
    )
    failed = samples.__class__.from_records(
        records=[*set(s.record for s in samples.failed)], config=config
    )

    if not complete and not failed:
        logger.info("No samples to update")

    if complete:
        logger.info(f"Marking {len(complete)} samples as complete")
        for sample in complete:
            logger.debug(f"Marking {sample.id} as complete")
        complete.set_state("complete")
        complete.update_derived(config)

    if failed:
        logger.warning(f"Marking {len(failed)} samples as failed")
        for sample in failed:
            logger.debug(f"Marking {sample.id} as failed")
        failed.set_state("error")
        failed.update_derived(config)
