from __future__ import annotations

import argparse
import sys

from . import __version__
from .workflow import run_workflow


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="pops",
        description="Reference-independent, group-level prioritisation of metagenomic contigs.",
    )
    p.add_argument("--group", required=True, help="Group sample manifest TSV")
    p.add_argument("--control", required=True, help="Control sample manifest TSV")
    p.add_argument("-x", type=int, default=2, help="Minimum number of group samples detecting a contig (default: 2)")
    p.add_argument("-y", type=int, default=1, help="Maximum number of control samples detecting a contig (default: 1)")
    p.add_argument("-outprefix", required=True, help="Output prefix")
    p.add_argument("-cfg", required=True, help="YAML configuration file")
    p.add_argument("--resume", action="store_true", help="Resume an interrupted identical run")
    p.add_argument("--keep-bam", action="store_true", help="Keep per-sample BAM and BAI files")
    p.add_argument("--version", action="version", version=f"POPS {__version__}")
    return p


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        run_workflow(
            group_manifest=args.group,
            control_manifest=args.control,
            x=args.x,
            y=args.y,
            outprefix=args.outprefix,
            cfg_path=args.cfg,
            resume=args.resume,
            keep_bam=args.keep_bam,
        )
    except Exception as exc:
        print(f"POPS ERROR: {exc}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
