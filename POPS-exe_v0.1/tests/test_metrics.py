from pops.metrics import interval_overlap_fraction, read_identity, reference_span_from_cigar


def test_interval_overlap_fraction():
    assert interval_overlap_fraction((0, 100), (80, 180)) == 0.2
    assert interval_overlap_fraction((0, 100), (20, 80)) == 1.0
    assert interval_overlap_fraction((0, 100), (100, 200)) == 0.0


def test_reference_span_from_cigar():
    assert reference_span_from_cigar("50M") == 50
    assert reference_span_from_cigar("10S40M") == 40
    assert reference_span_from_cigar("20M3D27M") == 50
    assert reference_span_from_cigar("10M100N40M") == 150


def test_read_identity_counts_paired_mates_separately():
    assert read_identity("frag001", 0x1 | 0x40) == ("frag001", 1, "")
    assert read_identity("frag001", 0x1 | 0x80) == ("frag001", 2, "")
    assert read_identity("frag001", 0, "") == ("frag001", 0, "")


def test_read_identity_secondary_alignment_same_mate_is_same_read():
    primary = read_identity("frag001", 0x1 | 0x40)
    secondary = read_identity("frag001", 0x1 | 0x40 | 0x100)
    assert primary == secondary


def test_read_identity_distinguishes_libraries():
    assert read_identity("frag001", 0, "DNA") != read_identity("frag001", 0, "RNA")
