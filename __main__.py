from cellophane import cellophane
from pathlib import Path


if __name__ == "__main__":
    _main = cellophane("QD-RNA", root=Path(__file__).parent)
    _main(prog_name="qd_rna")

