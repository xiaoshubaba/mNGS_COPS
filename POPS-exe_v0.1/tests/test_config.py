from pathlib import Path

from pops.config import load_config


def test_default_config_loads():
    path = Path(__file__).parents[1] / "config" / "default.yaml"
    cfg, raw = load_config(path)
    assert cfg.min_contig_length == 200
    assert cfg.mmseqs_min_seq_id == 0.90
    assert cfg.detection_max_overlap_fraction == 0.30
    assert raw["scientific"]["bowtie2_report_mode"] == "all"
