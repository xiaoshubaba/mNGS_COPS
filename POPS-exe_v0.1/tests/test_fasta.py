from pathlib import Path

from pops.fasta import fasta_lengths, filter_fasta_by_length, write_ranked_fasta


def test_length_filter_and_ranked_fasta(tmp_path: Path):
    src = tmp_path / "in.fa"
    src.write_text(">a\nAAAAA\n>b description\nAAA\n>c\nAAAAAAAA\n")
    dst = tmp_path / "filtered.fa"
    report = tmp_path / "report.tsv"
    n = filter_fasta_by_length(src, dst, report, 5)
    assert n == 2
    assert fasta_lengths(dst) == {"a": 5, "c": 8}
    ranked = tmp_path / "ranked.fa"
    write_ranked_fasta(dst, ranked, ["c", "a"])
    assert ranked.read_text().startswith(">c\nAAAAAAAA\n>a\nAAAAA\n")
