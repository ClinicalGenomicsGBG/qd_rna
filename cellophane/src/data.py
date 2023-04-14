"""Utilities for interacting with SLIMS"""

from collections import UserDict, UserList
from copy import deepcopy
from functools import reduce
from pathlib import Path
from typing import (
    Any,
    Callable,
    Hashable,
    Mapping,
    Optional,
    Sequence,
    Literal,
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
            case *k,:
                reduce(lambda d, k: d.setdefault(k, Container()), k[:-1], self.data)[
                    k[-1]
                ] = item
            case _:
                raise TypeError("Key must be a string or a sequence of strings")

    def __getitem__(self, key: Hashable | Sequence[Hashable]) -> Any:
        match key:
            case *k,:
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
        super().__init__(id=id, files=files, done=done, **kwargs)


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

    def split(self, link_by: Optional[str]) -> "Samples":
        if link_by is not None:
            linked = {s[link_by]: [l for l in self if l[link_by] == s[link_by]] for s in self}
            for linked in linked.values():
                yield self.__class__(linked)
        else:
            for sample in self:
                yield self.__class__([sample])

    def validate(self):
        for sample in self:
            if None in sample.files or not isinstance(sample.id, str):
                yield sample
        self.data = [s for s in self if None not in s.files]

    def __reduce__(self) -> Callable | tuple:
        return self.__class__, (self.data,)


class Mixin:
    sample_mixin: Optional[type]

    """Mixin class for adding properties to Samples"""
    def __init_subclass__(cls, sample_mixin: Optional[type] = None) -> None:
        cls.sample_mixin = sample_mixin
