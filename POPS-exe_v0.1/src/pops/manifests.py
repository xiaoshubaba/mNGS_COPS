from __future__ import annotations

import csv
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable


class ManifestError(ValueError):
    pass


@dataclass(frozen=True)
class Library:
    library_id: str
    read1: Path
    read2: Path | None
    metadata: dict[str, str]


@dataclass
class Sample:
    sample_id: str
    libraries: list[Library] = field(default_factory=list)

    @property
    def read1(self) -> Path:
        """Backward-compatible accessor for single-library samples."""
        if len(self.libraries) != 1:
            raise AttributeError("sample has multiple libraries; use sample.libraries")
        return self.libraries[0].read1

    @property
    def read2(self) -> Path | None:
        """Backward-compatible accessor for single-library samples."""
        if len(self.libraries) != 1:
            raise AttributeError("sample has multiple libraries; use sample.libraries")
        return self.libraries[0].read2

    @property
    def metadata(self) -> dict[str, str]:
        """Backward-compatible accessor for single-library samples."""
        if len(self.libraries) != 1:
            raise AttributeError("sample has multiple libraries; metadata are library-level")
        return self.libraries[0].metadata


def _resolve_read(value: str, manifest_dir: Path, required: bool) -> Path | None:
    value = (value or "").strip()
    if value in {"", "."}:
        if required:
            raise ManifestError("read1 is required for every library")
        return None
    p = Path(value).expanduser()
    if not p.is_absolute():
        p = (manifest_dir / p).resolve()
    else:
        p = p.resolve()
    if not p.is_file():
        raise ManifestError(f"FASTQ not found: {p}")
    return p


def read_manifest(path: str | Path) -> tuple[list[Sample], list[str]]:
    """Read a sample/library manifest.

    One row represents one sequencing library.  `sample_id` may repeat so that a
    biological sample can contribute multiple libraries (for example DNA and RNA).
    If `library_id` is omitted, each sample must occupy exactly one row and the
    library id defaults to the sample id for backward compatibility.
    """
    path = Path(path).resolve()
    if not path.is_file():
        raise ManifestError(f"Manifest not found: {path}")
    with path.open(newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        if reader.fieldnames is None:
            raise ManifestError(f"Manifest has no header: {path}")
        required = {"sample_id", "read1", "read2"}
        missing = required - set(reader.fieldnames)
        if missing:
            raise ManifestError(f"Manifest missing required columns: {', '.join(sorted(missing))}")
        has_library_id = "library_id" in reader.fieldnames
        sample_map: dict[str, Sample] = {}
        sample_order: list[str] = []
        seen_library_ids: set[str] = set()
        row_count_by_sample: dict[str, int] = {}

        for lineno, row in enumerate(reader, start=2):
            sid = (row.get("sample_id") or "").strip()
            if not sid:
                raise ManifestError(f"Empty sample_id at {path}:{lineno}")
            row_count_by_sample[sid] = row_count_by_sample.get(sid, 0) + 1
            if not has_library_id and row_count_by_sample[sid] > 1:
                raise ManifestError(
                    f"sample_id {sid!r} occurs on multiple rows in {path}; "
                    "add a unique library_id column when a biological sample has multiple libraries"
                )

            if has_library_id:
                lid = (row.get("library_id") or "").strip()
                if not lid:
                    raise ManifestError(f"Empty library_id at {path}:{lineno}")
            else:
                lid = sid
            if lid in seen_library_ids:
                raise ManifestError(f"Duplicate library_id {lid!r} in {path}")
            seen_library_ids.add(lid)

            r1 = _resolve_read(row.get("read1", ""), path.parent, True)
            r2 = _resolve_read(row.get("read2", ""), path.parent, False)
            structural = required | ({"library_id"} if has_library_id else set())
            metadata = {k: (v or "") for k, v in row.items() if k not in structural}
            lib = Library(lid, r1, r2, metadata)
            if sid not in sample_map:
                sample_map[sid] = Sample(sid)
                sample_order.append(sid)
            sample_map[sid].libraries.append(lib)

    if not sample_order:
        raise ManifestError(f"Manifest contains no samples: {path}")
    return [sample_map[sid] for sid in sample_order], list(reader.fieldnames)


def validate_unique_across(group: Iterable[Sample], control: Iterable[Sample]) -> None:
    g = {s.sample_id for s in group}
    c = {s.sample_id for s in control}
    overlap = sorted(g & c)
    if overlap:
        raise ManifestError("sample_id values must be unique across manifests; duplicated: " + ", ".join(overlap))

    glibs = {lib.library_id for s in group for lib in s.libraries}
    clibs = {lib.library_id for s in control for lib in s.libraries}
    loverlap = sorted(glibs & clibs)
    if loverlap:
        raise ManifestError("library_id values must be unique across manifests; duplicated: " + ", ".join(loverlap))
