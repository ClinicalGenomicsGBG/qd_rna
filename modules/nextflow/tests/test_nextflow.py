from pathlib import Path

from cellophane.src.testing import parametrize_from_yaml

ROOT = Path(__file__).parent

class Test_integration:
    @staticmethod
    @parametrize_from_yaml([ROOT / "integration.yaml"])
    def test_integration(definition: Path, run_definition):
        run_definition(definition)
