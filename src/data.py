"""Utilities for interacting with SLIMS"""

from collections import UserDict, UserList
from functools import reduce
from pathlib import Path
from typing import (
    Any,
    Hashable,
    Mapping,
    Optional,
    Sequence,
    TypeVar,
)

from yaml import safe_load


class Container(UserDict):
    """A dict that allows attribute access to its items"""

    def __contains__(self, key: Hashable | Sequence[Hashable]) -> bool:
        try:
            self[key]
        except (KeyError, TypeError):
            return False
        else:
            return True

    def __setitem__(self, key: Hashable | Sequence[Hashable], item: Any) -> None:
        if isinstance(item, Mapping) and not isinstance(item, Container):
            item = Container(item)

        match key:
            case k if isinstance(k, Hashable):
                self.data[k] = item
            case *k, :
                reduce(lambda d, k: d.setdefault(k, Container()), k[:-1], self.data)[
                    k[-1]
                ] = item
            case _:
                raise TypeError("Key must be a string or a sequence of strings")

    def __getitem__(self, key: Hashable | Sequence[Hashable]) -> Any:
        match key:
            case *k, :
                return reduce(lambda d, k: d[k], k, self.data)
            case k if isinstance(k, Hashable):
                return self.data[k]
            case k:
                raise TypeError("Key {k} is not hashble or a sequence of hashables")

    def __getattr__(self, key: str) -> Any:
        if key in dir(self):
            super().__getattribute__(key)
        elif "data" in self.__dict__ and key in self.data:
            return self.data[key]
        else:
            raise AttributeError(
                f"'{self.__class__.__name__}' object has no attribute '{key}'"
            )

    def __setattr__(self, key: str, value: Any) -> None:
        if key in dir(self) or key == "data":
            super().__setattr__(key, value)
        else:
            self[key] = value


class Sample(Container):
    """A basic sample container"""

    id: str
    files: Optional[list[str]]
    done: Optional[bool]

    def __init__(self, /, id, files=None, done=None, **kwargs):
        super().__init__(id=str(id), files=files, done=done, **kwargs)

    def __str__(self):
        return self.id


S = TypeVar("S", bound=Sample)


class Samples(UserList[S]):
    """A list of sample containers"""

    sample_class: type = Sample

    @classmethod
    def from_file(cls, path: Path):
        """Get samples from a YAML file"""
        with open(path, "r", encoding="utf-8") as handle:
            samples = []
            for sample in safe_load(handle):
                id = sample.pop("id")
                samples.append(Sample(id=str(id), **sample))
        return cls(samples)

    def split(self, link_by: Optional[str]):
        if link_by is not None:
            linked = {
                sample[link_by]: [
                    li for li in self if li[link_by] == sample[link_by]
                ]
                for sample in self
            }
            for li in linked.values():
                yield self.__class__(li)
        else:
            for sample in self:
                yield self.__class__([sample])

    def validate(self):
        for sample in self:
            if (
                sample.files is None
                or None in sample.files
                or not isinstance(sample.id, str)
            ):
                yield sample
        self.data = [
            sample
            for sample in self
            if sample.files is not None and None not in sample.files
        ]

    @property
    def unique_ids(self):
        return set(s.id for s in self)

    @property
    def complete(self):
        return self.__class__(
            [
                sample
                for samples_by_id in self.split(link_by="id")
                if all(sample.done for sample in samples_by_id)
                for sample in samples_by_id
            ]
        )

    @property
    def failed(self):
        return self.__class__(
            [
                sample
                for samples_by_id in self.split(link_by="id")
                if not all(sample.done for sample in samples_by_id)
                for sample in samples_by_id
            ]
        )

    def __reduce__(self):
        return self.__class__, (self.data,)

    def __str__(self):
        return "\n".join([str(s) for s in self])


class Mixin:
    sample_mixin: Optional[type]

    """Mixin class for adding properties to Samples"""

    def __init_subclass__(cls, sample_mixin: Optional[type] = None) -> None:
        cls.sample_mixin = sample_mixin
