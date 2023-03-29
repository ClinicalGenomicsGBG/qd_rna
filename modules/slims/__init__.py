"""Module for getting samples from SLIMS"""

from copy import deepcopy
from functools import cached_property, reduce
from json import loads
from logging import LoggerAdapter
from time import time
from typing import Any
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


def _parse_criteria(criteria: str | list[str]) -> list[Criterion]:
    """Parse criteria"""

    match criteria:
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
    *args: Criterion,
    connection: Slims,
    slims_id: str | list[str] | None = None,
    content_type: int | list[int] | None = None,
    max_age: int | str | None = None,
    derived_from: str | list[str] | None = None,
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

    match content_type:
        case int():
            criteria = criteria.add(equals("cntn_fk_contentType", content_type))
        case [*types]:
            criteria = criteria.add(is_one_of("cntn_fk_contentType", types))
        case _ if content_type is not None:
            raise TypeError(f"Invalid type for content_type: {type(content_type)}")

    match max_age:
        case int() | str():
            min_mtime = int(time() - parse_timespan(str(max_age))) * 1e3
            criteria = criteria.add(greater_than("cntn_modifiedOn", min_mtime))
        case _ if max_age is not None:
            raise TypeError(f"Expected int or str, got {type(max_age)}")

    match derived_from:
        case None:
            pass
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

    record: Record | None
    bioinformatics: Record | None

    @classmethod
    def from_record(cls, record: Record, config: cfg.Config, **kwargs):
        """Create a sample from a SLIMS fastq record"""
        return cls(
            id=record.cntn_id.value,
            bioinformatics=None,
            record=record,
            **{
                key: _get_field(record, field)
                for key, field in config.slims.map_field.items()
            },
            **kwargs,
        )

    def add_bioinformatics(
        self,
        config: cfg.Config,
    ):
        """Add a bioinformatics record to the sample"""

        if (
            "bioinformatics" in self
            and self.bioinformatics is None
            and self._connection is not None
        ):
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
                if "bioinformatics" in self and self.bioinformatics is not None:
                    self.bioinformatics = self.bioinformatics.update(
                        {config.slims.bioinfo.state_field: state}
                    )
            case _:
                raise ValueError(f"Invalid state: {state}")

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
) -> SlimsSamples | None:
    """Load novel samples from SLIMS."""
    if "slims" in config:
        slims_connection = Slims(
            name=__package__,
            url=config.slims.url,
            username=config.slims.username,
            password=config.slims.password,
        )

        if samples:
            logger.info("Augmenting existing samples with info from SLIMS")
            slims_ids = [s.id for s in samples]
            max_age = None

        elif "id" in config.slims:
            logger.info("Fetching samples from SLIMS by ID")
            slims_ids = config.slims.id
            max_age = None

        else:
            logger.info(f"Fetching samples from the last {config.slims.novel_max_age}")
            slims_ids = None
            max_age: str = config.slims.novel_max_age

        criteria = _parse_criteria(config.slims.criteria)

        parent_records: list[Record] | None = None
        for criterion in criteria[:-1]:
            parent_records = get_records(
                criterion,
                connection=slims_connection,
                derived_from=parent_records,
                slims_id=slims_ids if not config.slims.unrestrict_parents else None,
                max_age=max_age if not config.slims.unrestrict_parents else None,
            )

        if parent_records:
            logger.debug(
                f"Found parent records: {[r.cntn_id.value for r in parent_records]}"
            )

        records = get_records(
            criteria[-1],
            connection=slims_connection,
            derived_from=parent_records,
            slims_id=slims_ids,
            max_age=max_age,
        )

        if config.slims.bioinfo.check:
            logger.info("Checking SLIMS for completed bioinformatics")
            bioinfo = get_records(
                _parse_criteria(config.slims.bioinfo.check_criteria),
                connection=slims_connection,
                derived_from=records,
                content_type=config.slims.bioinfo.content_type,
            )

            original_ids = [r.cntn_id.value for r in records]
            records = [
                record
                for record in records
                if record.pk() not in [b.cntn_fk_originalContent.value for b in bioinfo]
            ]

            for sid in set(original_ids) - set([r.cntn_id.value for r in records]):
                logger.info(f"Found completed bioinformatics for {sid}")

        slims_samples = samples.from_records(records, config)

        if samples:
            for idx, sample in enumerate(samples):
                match = [m for m in slims_samples if m.id == sample.id]
                common_keys = set([k for s in match for k in s]) & set(sample.keys())
                common_keys -= set(["files", "backup", "complete"])
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
    """Add bioinformatics content to SLIMS samples"""
    if config.slims.dry_run:
        logger.debug("Dry run - Not adding bioinformatics")
    elif config.slims.bioinfo.create and samples:
        logger.info("Creating bioinformatics records")
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
        unique = {
            pk: [s for s in samples if s.record.pk() == pk]
            for pk in set(s.record.pk() for s in samples)
        }
        for _samples in unique.values():
            if all(s.complete for s in _samples):
                logger.info(f"Marking {len(unique)} samples as complete")
                _samples[0].set_bioinformatics_state("complete", config)
            else:
                logger.warning(f"Marking {len(unique)} samples as failed")
                _samples[0].set_bioinformatics_state("error", config)
    else:
        logger.info("No SLIMS samples to update")
