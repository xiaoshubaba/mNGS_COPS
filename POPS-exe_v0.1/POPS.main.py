#!/usr/bin/env python3
"""Compatibility entry point for running POPS directly from the source tree."""
from __future__ import annotations

import sys
from pathlib import Path

# Support the documented source-tree invocation:
#   python3 POPS.main.py ...
# without requiring an editable/package installation first.
_REPO_ROOT = Path(__file__).resolve().parent
_SRC = _REPO_ROOT / "src"
if _SRC.is_dir():
    sys.path.insert(0, str(_SRC))

from pops.cli import main


if __name__ == "__main__":
    raise SystemExit(main())
