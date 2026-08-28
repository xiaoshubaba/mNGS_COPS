from __future__ import annotations

import subprocess
import sys
from pathlib import Path


def test_source_tree_entrypoint_help_runs_without_install():
    root = Path(__file__).resolve().parents[1]
    proc = subprocess.run(
        [sys.executable, str(root / "POPS.main.py"), "--help"],
        cwd=root,
        text=True,
        capture_output=True,
        check=False,
    )
    assert proc.returncode == 0, proc.stderr
    assert "--group" in proc.stdout
    assert "--control" in proc.stdout


def test_cli_x_y_defaults():
    from pops.cli import build_parser
    parser = build_parser()
    args = parser.parse_args([
        "--group", "g.tsv", "--control", "c.tsv",
        "-outprefix", "out/run", "-cfg", "config/default.yaml"
    ])
    assert args.x == 2
    assert args.y == 1
