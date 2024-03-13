from pathlib import Path
from unittest.mock import MagicMock

from cellophane import data
from cellophane.src.testing import parametrize_from_yaml
from pytest import mark, param, raises
from ruamel.yaml import YAML
from slims.slims import Record
from slims_.src import util as slims_util


class Test_integration:
    @staticmethod
    @parametrize_from_yaml("integration.yaml")
    def test_integration(definition: Path, run_definition):
        run_definition(definition)


class RecordMock(MagicMock):
    def __init__(self, **kwargs):
        super().__init__(spec=Record)
        self.__dict__.update(
            {
                k: data.Container(v) if isinstance(v, dict) else v
                for k, v in kwargs.items()
            }
        )
    def __reduce__(self):
        return (RecordMock, ())

class Test_criteria:
    @staticmethod
    @mark.parametrize(
        "criteria,slims,exception,kwargs",
        [
            param(
                d["criteria"],
                d.get("slims"),
                d.get("exception"),
                d.get("kwargs", {}),
                id=d["id"],
            )
            for d in YAML(typ="unsafe").load_all("criteria.yaml")
        ],
    )
    def test_criteria(criteria, slims, exception, kwargs):
        if exception:
            with raises(exception):
                slims_util.parse_criteria(criteria, **kwargs)
        else:
            parsed = slims_util.parse_criteria(criteria, **kwargs)
            assert [c.to_dict() for c in parsed] == slims
