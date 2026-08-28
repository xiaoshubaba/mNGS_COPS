from pathlib import Path

import pytest

from pops.manifests import ManifestError, read_manifest, validate_unique_across


def _fq(path: Path):
    path.write_text("@r1\nA\n+\nI\n")


def test_manifest_resolves_relative_paths_and_metadata(tmp_path: Path):
    _fq(tmp_path / "a.fq")
    manifest = tmp_path / "group.tsv"
    manifest.write_text("sample_id\tread1\tread2\tnote\nS1\ta.fq\t.\tx\n")
    samples, fields = read_manifest(manifest)
    assert samples[0].read1 == (tmp_path / "a.fq").resolve()
    assert samples[0].read2 is None
    assert samples[0].metadata == {"note": "x"}
    assert fields == ["sample_id", "read1", "read2", "note"]


def test_unique_across_manifests(tmp_path: Path):
    _fq(tmp_path / "a.fq")
    g = tmp_path / "g.tsv"
    c = tmp_path / "c.tsv"
    g.write_text("sample_id\tread1\tread2\nS1\ta.fq\t.\n")
    c.write_text("sample_id\tread1\tread2\nS1\ta.fq\t.\n")
    gs, _ = read_manifest(g)
    cs, _ = read_manifest(c)
    with pytest.raises(ManifestError):
        validate_unique_across(gs, cs)


def test_manifest_allows_multiple_libraries_per_sample_with_library_id(tmp_path: Path):
    _fq(tmp_path / "dna.fq")
    _fq(tmp_path / "rna.fq")
    manifest = tmp_path / "group.tsv"
    manifest.write_text(
        "sample_id\tlibrary_id\tread1\tread2\tlibrary_type\n"
        "S1\tS1_DNA\tdna.fq\t.\tDNA\n"
        "S1\tS1_RNA\trna.fq\t.\tRNA\n"
    )
    samples, fields = read_manifest(manifest)
    assert len(samples) == 1
    assert samples[0].sample_id == "S1"
    assert [x.library_id for x in samples[0].libraries] == ["S1_DNA", "S1_RNA"]
    assert [x.metadata["library_type"] for x in samples[0].libraries] == ["DNA", "RNA"]
    assert fields == ["sample_id", "library_id", "read1", "read2", "library_type"]


def test_repeated_sample_requires_library_id(tmp_path: Path):
    _fq(tmp_path / "a.fq")
    _fq(tmp_path / "b.fq")
    manifest = tmp_path / "group.tsv"
    manifest.write_text(
        "sample_id\tread1\tread2\n"
        "S1\ta.fq\t.\n"
        "S1\tb.fq\t.\n"
    )
    with pytest.raises(ManifestError):
        read_manifest(manifest)
