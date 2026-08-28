from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path


class RankingError(RuntimeError):
    pass


@dataclass
class RankedContig:
    contig_id: str
    length: int
    group_detected_n: int
    control_detected_n: int
    group_mean_depth: float
    control_mean_depth: float
    denominator_depth: float
    score: float
    raw_rank: int = 0
    pass_x: bool = False
    pass_y: bool = False
    retained_rank: int | None = None
    exclusion_reason: str = ""


def read_numeric_matrix(path: str | Path) -> tuple[list[str], dict[str, list[float]]]:
    with Path(path).open(newline="") as fh:
        reader = csv.reader(fh, delimiter="\t")
        header = next(reader)
        samples = header[1:]
        data = {row[0]: [float(x) for x in row[1:]] for row in reader}
    return samples, data


def calculate_rankings(
    contig_lengths: dict[str, int],
    group_depth: dict[str, list[float]],
    control_depth: dict[str, list[float]],
    group_detected: dict[str, list[float]],
    control_detected: dict[str, list[float]],
    x: int,
    y: int,
) -> list[RankedContig]:
    control_means = {
        cid: (sum(control_depth[cid]) / len(control_depth[cid])) for cid in contig_lengths
    }
    positive = [v for v in control_means.values() if v > 0]
    if not positive:
        raise RankingError(
            "All representative contigs have zero mean depth in all controls; "
            "the prespecified zero-denominator replacement is undefined."
        )
    replacement = min(positive)
    rows: list[RankedContig] = []
    for cid, length in contig_lengths.items():
        gd = sum(group_depth[cid]) / len(group_depth[cid])
        cd = control_means[cid]
        denom = cd if cd > 0 else replacement
        gn = int(sum(group_detected[cid]))
        cn = int(sum(control_detected[cid]))
        px = gn >= x
        py = cn <= y
        reasons = []
        if not px:
            reasons.append("case_recurrence")
        if not py:
            reasons.append("control_frequency")
        rows.append(RankedContig(
            contig_id=cid,
            length=length,
            group_detected_n=gn,
            control_detected_n=cn,
            group_mean_depth=gd,
            control_mean_depth=cd,
            denominator_depth=denom,
            score=gd / denom,
            pass_x=px,
            pass_y=py,
            exclusion_reason=";".join(reasons),
        ))

    rows.sort(key=lambda r: (-r.score, -r.group_mean_depth, r.contig_id))
    for i, row in enumerate(rows, 1):
        row.raw_rank = i
    retained = [r for r in rows if r.pass_x and r.pass_y]
    retained.sort(key=lambda r: (-r.score, -r.group_mean_depth, r.contig_id))
    for i, row in enumerate(retained, 1):
        row.retained_rank = i
    return rows


def write_ranking_tables(rows: list[RankedContig], all_path: str | Path, ranked_path: str | Path) -> list[str]:
    cols = [
        "contig_id", "length", "raw_rank", "retained_rank", "pass_x", "pass_y",
        "group_detected_n", "control_detected_n", "group_mean_depth",
        "control_mean_depth", "denominator_depth", "score", "exclusion_reason",
    ]
    with Path(all_path).open("w", newline="") as out:
        w = csv.writer(out, delimiter="\t")
        w.writerow(cols)
        for r in sorted(rows, key=lambda z: z.raw_rank):
            w.writerow([
                r.contig_id, r.length, r.raw_rank, r.retained_rank or "", int(r.pass_x), int(r.pass_y),
                r.group_detected_n, r.control_detected_n,
                f"{r.group_mean_depth:.12g}", f"{r.control_mean_depth:.12g}",
                f"{r.denominator_depth:.12g}", f"{r.score:.12g}", r.exclusion_reason,
            ])
    retained = sorted((r for r in rows if r.retained_rank is not None), key=lambda z: z.retained_rank)
    with Path(ranked_path).open("w", newline="") as out:
        w = csv.writer(out, delimiter="\t")
        w.writerow(cols)
        for r in retained:
            w.writerow([
                r.contig_id, r.length, r.raw_rank, r.retained_rank, int(r.pass_x), int(r.pass_y),
                r.group_detected_n, r.control_detected_n,
                f"{r.group_mean_depth:.12g}", f"{r.control_mean_depth:.12g}",
                f"{r.denominator_depth:.12g}", f"{r.score:.12g}", r.exclusion_reason,
            ])
    return [r.contig_id for r in retained]
