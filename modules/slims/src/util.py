"""Module for getting samples from SLIMS"""

import re
from datetime import datetime, timedelta
from functools import cache, reduce
from json import loads
from typing import Any

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


@cache
def split_criteria(criteria: str) -> list[str]:
    """
    Split string on "and"/"or" but not within parentheses

    >>> _parse_bool("a is x and (b is y or c is d) or g is h")
    ['a is x', 'and', 'b is y or c is d', 'or', 'g is h']
    """

    parts = []
    _criteria = " ".join(criteria.split())
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
        elif _criteria[0:5] == " and " and depth == 0:
            parts = [part, "and", _criteria[5:]]
            break
        elif _criteria[0:4] == " or " and depth == 0:
            parts = [part, "or", _criteria[4:]]
            break
        else:
            part += _criteria[0]

        _criteria = _criteria[1:]

    if depth != 0:
        raise ValueError(f"Unmatched parentheses: {criteria}")

    return parts or [part]


def parse_criteria(  # type: ignore[return]
    criteria: str | list[str],
    parent_records: list[Record] | None = None,
) -> list[Criterion]:
    """Parse criteria"""

    match criteria:
        case str(criteria) if all(w not in criteria for w in [
            "equals",
            "not_equals",
            "one_of",
            "not_one_of",
            "equals_ignore_case",
            "not_equals_ignore_case",
            "contains",
            "not_contains",
            "starts_with",
            "not_starts_with",
            "ends_with",
            "not_ends_with",
            "between",
            "not_between",
            "greater_than",
            "less_than",
        ]):
            raise ValueError(f"Invalid criteria: {criteria}")
        case str(criteria) if criteria.startswith("->"):
            if parent_records is None:
                raise ValueError("Cannot use leading '->' without parent record(s)")
            else:
                _parsed = parse_criteria(criteria[2:])
                parent_pks = [p.pk() for p in parent_records]
                return [
                    conjunction()
                    .add(is_one_of("cntn_fk_originalContent", parent_pks))
                    .add(_parsed[0]),
                    *_parsed[1:],
                ]
        case str(criteria) if parent_records:
            _parsed = parse_criteria(criteria)
            parent_pks = [p.pk() for p in parent_records]
            return [
                conjunction().add(is_one_of("cntn_pk", parent_pks)).add(_parsed[0]),
                *_parsed[1:],
            ]
        case str(criteria) if "->" in criteria:
            return [parse_criteria(c)[0] for c in criteria.split("->")]
        case str(criteria):
            return parse_criteria(split_criteria(criteria))

        case [criterion]:
            criteria = split_criteria(criterion)
            if len(criteria) == 1:
                return parse_criteria(criteria[0].split(" "))
            else:
                return parse_criteria(criteria)

        case [a, "and", b]:
            return [conjunction().add(parse_criteria(a)[0]).add(parse_criteria(b)[0])]
        case [a, "or", b]:
            return [disjunction().add(parse_criteria(a)[0]).add(parse_criteria(b)[0])]

        case [_, *mid, _] if "and" in mid or "or" in mid:
            # This handles cases where multiple parentheses are used
            # around a single criterion (e.g. "((((a equals b))))")
            return parse_criteria(" ".join(criteria))

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
        case [field, "not_equals_ignore_case", value]:
            return [is_not(equals_ignore_case(field, value))]

        case [field, "contains", value]:
            return [contains(field, value)]
        case [field, "not_contains", value]:
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


def get_field(record: Record, field: str, default=None) -> Any:
    """Get a field from SLIMS record"""
    try:
        if field.startswith("json:"):
            _field, *_key = re.split(r"\.|(\[[0-9]*\])", field[5:])
            _key = [int(k.strip("[]")) if k.startswith("[") else k for k in _key if k]
            _json = loads(record.__dict__[_field].value)
            return reduce(lambda x, y: x[y], _key, _json)
        else:
            return getattr(record, field).value
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
        _parsed = parse_criteria(string_criteria, parent_records=derived_from)
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

