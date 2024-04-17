from pathlib import Path
from unittest.mock import MagicMock

from cellophane import Sample
from cellophane.src.testing import parametrize_from_yaml
from pytest import fixture, mark, param, raises
from pytest_mock import MockerFixture

from .. import HCPSample, fetch

ROOT = Path(__file__).parent

class Test_integration:
    @staticmethod
    @parametrize_from_yaml([ROOT / "integration.yaml"])
    def test_integration(definition: Path, run_definition):
        run_definition(definition)

