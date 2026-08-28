# Contributing

Changes that alter scientific behaviour should be explicit, tested, documented,
and made through a new configuration/output prefix when reproducing manuscript
analyses. Do not add taxonomic information to filtering, scoring, or ranking.

Before opening a pull request:

```bash
python -m pip install -e . pytest
pytest
```

Changes to external-tool command lines should also update `USAGE.md` and, when
relevant, `config/default.yaml`.
