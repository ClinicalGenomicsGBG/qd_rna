import sys
from cellophane import cellophane
from pathlib import Path


if __name__ == "__main__":
    _root = Path(__file__).parent
    sys.path.append(str(_root))
    _main = cellophane("QD-RNA", root=_root)
    _main(prog_name="qd_rna")
