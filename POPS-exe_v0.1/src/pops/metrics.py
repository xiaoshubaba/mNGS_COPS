from __future__ import annotations

import csv
import re
import subprocess
from pathlib import Path

from .fasta import fasta_lengths

_CIGAR_RE = re.compile(r"(\d+)([MIDNSHP=X])")


def interval_overlap_fraction(a: tuple[int, int], b: tuple[int, int]) -> float:
    a0, a1 = a
    b0, b1 = b
    la = max(0, a1 - a0)
    lb = max(0, b1 - b0)
    shorter = min(la, lb)
    if shorter == 0:
        return 1.0
    overlap = max(0, min(a1, b1) - max(a0, b0))
    return overlap / shorter


def reference_span_from_cigar(cigar: str) -> int:
    """Reference-consuming alignment span for a SAM CIGAR string."""
    if cigar == "*":
        return 0
    span = 0
    consumed = 0
    for n_s, op in _CIGAR_RE.findall(cigar):
        n = int(n_s)
        consumed += len(n_s) + 1
        if op in {"M", "D", "N", "=", "X"}:
            span += n
    if consumed != len(cigar):
        raise ValueError(f"Invalid CIGAR string: {cigar}")
    return span


def read_identity(qname: str, flag: int, read_group: str = "") -> tuple[str, int, str]:
    """Return a stable identity for one sequenced read.

    Paired-end mates share a QNAME in SAM/BAM, but R1 and R2 are distinct reads
    for the POPS detection rule. The second tuple element is 1 for first mate,
    2 for second mate, and 0 for an unpaired/unspecified read. The read-group
    element distinguishes identical QNAMEs arising from different sequencing
    libraries of the same biological sample. Multiple alignments of the same
    mate within one library therefore retain the same identity.
    """
    if flag & 0x40:
        mate = 1
    elif flag & 0x80:
        mate = 2
    else:
        mate = 0
    return qname, mate, read_group


def has_two_distinct_reads(bam: str | Path, contig: str, max_overlap_fraction: float) -> bool:
    """Test the two-distinct-read rule from `samtools view BAM contig` output.

    Paired-end R1 and R2 count as two distinct reads even when they share the
    same QNAME. Multiple alignments from the same mate cannot by themselves
    satisfy the criterion.
    """
    proc = subprocess.Popen(
        ["samtools", "view", str(bam), contig],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )
    assert proc.stdout is not None
    intervals: list[tuple[tuple[str, int, str], tuple[int, int]]] = []
    seen_reads: set[tuple[str, int, str]] = set()
    try:
        for line in proc.stdout:
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 6:
                continue
            qname = cols[0]
            flag = int(cols[1])
            read_group = ""
            for tag in cols[11:]:
                if tag.startswith("RG:Z:"):
                    read_group = tag[5:]
                    break
            rid = read_identity(qname, flag, read_group)
            if flag & 0x4 or rid in seen_reads:
                continue
            pos0 = int(cols[3]) - 1
            span = reference_span_from_cigar(cols[5])
            if span <= 0:
                continue
            interval = (pos0, pos0 + span)
            for prev_rid, prev_interval in intervals:
                if prev_rid != rid and interval_overlap_fraction(interval, prev_interval) < max_overlap_fraction:
                    proc.stdout.close()
                    proc.terminate()
                    proc.wait(timeout=5)
                    return True
            intervals.append((rid, interval))
            seen_reads.add(rid)
    finally:
        if proc.stdout and not proc.stdout.closed:
            proc.stdout.close()
    stderr = proc.stderr.read() if proc.stderr else ""
    rc = proc.wait()
    if rc != 0:
        raise RuntimeError(f"samtools view failed for {contig}: {stderr[-2000:]}")
    return False


def parse_bedtools_coverage(path: str | Path) -> dict[str, float]:
    result: dict[str, float] = {}
    with Path(path).open() as fh:
        for line in fh:
            if not line.strip():
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 7:
                raise ValueError(f"Unexpected bedtools coverage output: {line.rstrip()}")
            result[cols[3]] = float(cols[-1])
    return result


def parse_bedtools_mean(path: str | Path) -> dict[str, float]:
    result: dict[str, float] = {}
    with Path(path).open() as fh:
        for line in fh:
            if not line.strip():
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 5:
                raise ValueError(f"Unexpected bedtools coverage -mean output: {line.rstrip()}")
            result[cols[3]] = float(cols[-1])
    return result


def write_contig_bed(representatives_fasta: str | Path, out_bed: str | Path) -> None:
    lengths = fasta_lengths(representatives_fasta)
    with Path(out_bed).open("w") as out:
        for contig, length in lengths.items():
            out.write(f"{contig}\t0\t{length}\t{contig}\n")


def write_sample_metrics(
    sample_id: str,
    representatives_fasta: str | Path,
    bam: str | Path,
    breadth_file: str | Path,
    depth_file: str | Path,
    out_tsv: str | Path,
    max_overlap_fraction: float,
) -> None:
    lengths = fasta_lengths(representatives_fasta)
    breadth = parse_bedtools_coverage(breadth_file)
    depth = parse_bedtools_mean(depth_file)
    with Path(out_tsv).open("w", newline="") as out:
        writer = csv.writer(out, delimiter="\t")
        writer.writerow(["contig_id", "sample_id", "breadth", "mean_depth", "detected"])
        for contig in lengths:
            detected = has_two_distinct_reads(bam, contig, max_overlap_fraction)
            writer.writerow([
                contig,
                sample_id,
                f"{breadth.get(contig, 0.0):.12g}",
                f"{depth.get(contig, 0.0):.12g}",
                1 if detected else 0,
            ])


def read_sample_metric(path: str | Path, field: str) -> dict[str, str]:
    with Path(path).open(newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        return {row["contig_id"]: row[field] for row in reader}


def merge_metric_matrix(
    contig_ids: list[str],
    sample_ids: list[str],
    sample_metric_paths: dict[str, Path],
    field: str,
    out_path: str | Path,
    batch_size: int = 128,
) -> None:
    """Merge per-sample metric files through temporary <=batch_size wide matrices."""
    out_path = Path(out_path)
    batch_paths: list[Path] = []
    for bi, start in enumerate(range(0, len(sample_ids), batch_size)):
        batch = sample_ids[start:start + batch_size]
        data = {sid: read_sample_metric(sample_metric_paths[sid], field) for sid in batch}
        bp = out_path.with_suffix(out_path.suffix + f".batch{bi:04d}.tmp")
        with bp.open("w") as out:
            out.write("contig_id\t" + "\t".join(batch) + "\n")
            for cid in contig_ids:
                out.write(cid + "\t" + "\t".join(data[sid].get(cid, "0") for sid in batch) + "\n")
        batch_paths.append(bp)

    handles = [p.open() for p in batch_paths]
    try:
        headers = [h.readline().rstrip("\n").split("\t") for h in handles]
        with out_path.open("w") as out:
            all_samples = [x for hdr in headers for x in hdr[1:]]
            out.write("contig_id\t" + "\t".join(all_samples) + "\n")
            for cid in contig_ids:
                rows = [h.readline().rstrip("\n").split("\t") for h in handles]
                if any(not row or row[0] != cid for row in rows):
                    raise ValueError(f"Metric merge lost contig order at {cid}")
                vals = [v for row in rows for v in row[1:]]
                out.write(cid + "\t" + "\t".join(vals) + "\n")
    finally:
        for h in handles:
            h.close()
        for p in batch_paths:
            p.unlink(missing_ok=True)
