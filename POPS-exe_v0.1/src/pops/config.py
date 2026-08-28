from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any

import yaml


class ConfigError(ValueError):
    pass


@dataclass(frozen=True)
class ToolConfig:
    min_contig_length: int
    mmseqs_min_seq_id: float
    mmseqs_min_coverage: float
    mmseqs_cov_mode: int
    bowtie2_preset: str
    bowtie2_report_mode: str
    detection_max_overlap_fraction: float
    threads_megahit: int
    threads_mmseqs: int
    threads_bowtie2: int
    mapping_jobs: int
    merge_batch_size: int


DEFAULTS = {
    "scientific": {
        "min_contig_length": 200,
        "mmseqs_min_seq_id": 0.90,
        "mmseqs_min_coverage": 0.90,
        "mmseqs_cov_mode": 0,
        "bowtie2_preset": "very-sensitive",
        "bowtie2_report_mode": "all",
        "detection_max_overlap_fraction": 0.30,
    },
    "resources": {
        "threads_megahit": 16,
        "threads_mmseqs": 8,
        "threads_bowtie2": 4,
        "mapping_jobs": 4,
        "merge_batch_size": 128,
    },
}


def _validate_fraction(name: str, value: Any) -> float:
    try:
        v = float(value)
    except (TypeError, ValueError) as exc:
        raise ConfigError(f"{name} must be numeric") from exc
    if not 0 <= v <= 1:
        raise ConfigError(f"{name} must be between 0 and 1")
    return v


def load_config(path: str | Path) -> tuple[ToolConfig, dict[str, Any]]:
    path = Path(path)
    if not path.is_file():
        raise ConfigError(f"Configuration file not found: {path}")
    raw = yaml.safe_load(path.read_text()) or {}
    if not isinstance(raw, dict):
        raise ConfigError("Configuration root must be a mapping")

    allowed_top = {"scientific", "resources"}
    unknown_top = set(raw) - allowed_top
    if unknown_top:
        raise ConfigError(f"Unknown configuration section(s): {', '.join(sorted(unknown_top))}")

    merged = {
        "scientific": {**DEFAULTS["scientific"], **(raw.get("scientific") or {})},
        "resources": {**DEFAULTS["resources"], **(raw.get("resources") or {})},
    }

    allowed_sci = set(DEFAULTS["scientific"])
    allowed_res = set(DEFAULTS["resources"])
    unknown_sci = set(merged["scientific"]) - allowed_sci
    unknown_res = set(merged["resources"]) - allowed_res
    if unknown_sci or unknown_res:
        bad = sorted(unknown_sci | unknown_res)
        raise ConfigError(f"Unknown configuration key(s): {', '.join(bad)}")

    sci = merged["scientific"]
    res = merged["resources"]
    if int(sci["min_contig_length"]) < 1:
        raise ConfigError("min_contig_length must be >=1")
    report_mode = str(sci["bowtie2_report_mode"]).lower()
    if report_mode != "all":
        raise ConfigError("Only bowtie2_report_mode: all is supported for manuscript-compatible runs")
    if int(res["merge_batch_size"]) < 1:
        raise ConfigError("merge_batch_size must be >=1")
    for key in ("threads_megahit", "threads_mmseqs", "threads_bowtie2", "mapping_jobs"):
        if int(res[key]) < 1:
            raise ConfigError(f"{key} must be >=1")

    cfg = ToolConfig(
        min_contig_length=int(sci["min_contig_length"]),
        mmseqs_min_seq_id=_validate_fraction("mmseqs_min_seq_id", sci["mmseqs_min_seq_id"]),
        mmseqs_min_coverage=_validate_fraction("mmseqs_min_coverage", sci["mmseqs_min_coverage"]),
        mmseqs_cov_mode=int(sci["mmseqs_cov_mode"]),
        bowtie2_preset=str(sci["bowtie2_preset"]),
        bowtie2_report_mode=report_mode,
        detection_max_overlap_fraction=_validate_fraction(
            "detection_max_overlap_fraction", sci["detection_max_overlap_fraction"]
        ),
        threads_megahit=int(res["threads_megahit"]),
        threads_mmseqs=int(res["threads_mmseqs"]),
        threads_bowtie2=int(res["threads_bowtie2"]),
        mapping_jobs=int(res["mapping_jobs"]),
        merge_batch_size=int(res["merge_batch_size"]),
    )
    return cfg, merged
