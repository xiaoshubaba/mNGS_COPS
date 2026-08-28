# Changelog

## 0.1.0

- Clarified and implemented the manuscript detection rule so paired-end R1 and R2 are counted as two distinct reads even when they share a QNAME; multiple alignments from the same mate are counted once.

- Initial public repository implementation of the workflow described in README.
- MEGAHIT group coassembly, length filtering, MMseqs2 representative clustering,
  Bowtie2 `-a` sample mapping, BEDTools breadth/depth metrics, two-distinct-read
  detection, X/Y filtering, ranking, provenance logging, Docker deployment, and
  interruption-safe resume support.
- CLI now defaults to `-x 2 -y 1`, matching the manuscript settings.
- Manifests now support multiple sequencing libraries per biological sample via repeated `sample_id` rows and `library_id`; X/Y recurrence is calculated at the biological-sample level.
- DNA and RNA libraries from the same specimen can be combined while retaining library-level provenance.

