from __future__ import annotations

from pathlib import Path
from typing import Iterator


def read_fasta(path: str | Path) -> Iterator[tuple[str, str]]:
    name = None
    chunks: list[str] = []
    with Path(path).open() as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(chunks)
                name = line[1:].split()[0]
                chunks = []
            else:
                if name is None:
                    raise ValueError(f"Invalid FASTA: sequence before header in {path}")
                chunks.append(line.strip())
    if name is not None:
        yield name, "".join(chunks)


def filter_fasta_by_length(source: str | Path, target: str | Path, report: str | Path, min_len: int) -> int:
    kept = 0
    with Path(target).open("w") as out, Path(report).open("w") as rep:
        rep.write("contig_id\tlength\tpasses_length_filter\n")
        for name, seq in read_fasta(source):
            ok = len(seq) >= min_len
            rep.write(f"{name}\t{len(seq)}\t{1 if ok else 0}\n")
            if ok:
                kept += 1
                out.write(f">{name}\n")
                for i in range(0, len(seq), 80):
                    out.write(seq[i:i+80] + "\n")
    return kept


def fasta_lengths(path: str | Path) -> dict[str, int]:
    return {name: len(seq) for name, seq in read_fasta(path)}


def write_ranked_fasta(source: str | Path, target: str | Path, ordered_ids: list[str]) -> None:
    seqs = dict(read_fasta(source))
    missing = [x for x in ordered_ids if x not in seqs]
    if missing:
        raise ValueError(f"Contig(s) missing from representative FASTA: {', '.join(missing[:5])}")
    with Path(target).open("w") as out:
        for name in ordered_ids:
            seq = seqs[name]
            out.write(f">{name}\n")
            for i in range(0, len(seq), 80):
                out.write(seq[i:i+80] + "\n")
