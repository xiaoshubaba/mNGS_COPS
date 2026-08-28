import pytest

from pops.ranking import RankingError, calculate_rankings


def test_ranking_zero_control_replacement_and_filters():
    lengths = {"a": 100, "b": 100, "c": 100}
    gdepth = {"a": [10, 10], "b": [5, 5], "c": [30, 30]}
    cdepth = {"a": [0, 0], "b": [1, 1], "c": [3, 3]}
    gdet = {"a": [1, 1], "b": [1, 0], "c": [1, 1]}
    cdet = {"a": [0, 0], "b": [0, 0], "c": [1, 1]}
    rows = calculate_rankings(lengths, gdepth, cdepth, gdet, cdet, x=2, y=1)
    by_id = {r.contig_id: r for r in rows}
    # Smallest non-zero control mean is 1, so a scores 10/1.
    assert by_id["a"].denominator_depth == 1
    assert by_id["a"].score == 10
    # c has two control detections and fails Y=1.
    assert not by_id["c"].pass_y
    # b fails recurrence X=2.
    assert not by_id["b"].pass_x
    assert by_id["a"].retained_rank == 1


def test_ranking_stops_when_all_control_depths_zero():
    with pytest.raises(RankingError):
        calculate_rankings(
            {"a": 100}, {"a": [1]}, {"a": [0]}, {"a": [1]}, {"a": [0]}, x=1, y=1
        )
