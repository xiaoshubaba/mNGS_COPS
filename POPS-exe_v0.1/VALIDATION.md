# Implementation validation notes

This repository implements the workflow defined in the project README and keeps
taxonomy outside filtering, scoring, and ranking.

## Implemented invariants

- Group FASTQs alone are supplied to MEGAHIT coassembly.
- Minimum contig length is applied before MMseqs2 clustering.
- MMseqs2 defaults are 0.90 identity, 0.90 coverage, and cov-mode 0.
- Every group and control sample is mapped independently to representative
  contigs with Bowtie2 `-a` and the configured preset.
- BEDTools reports per-sample breadth and mean depth.
- Binary detection requires two distinct reads with aligned-interval overlap
  below 0.30 of the shorter interval. Paired-end R1 and R2 count separately even
  when they share a QNAME; multiple alignments from the same mate count once.
- X and Y are applied to the binary detection matrices.
- Raw ranking precedes X/Y filtering; retained ranking follows filtering.
- Zero control mean depths use the smallest non-zero contig-level control mean
  depth from the same analysis; the workflow stops if that value does not exist.
- Ties are deterministic: score, then group mean depth, then contig identifier.
- Taxonomic information is never read by the workflow.
- A completed output prefix cannot be resumed or overwritten.

## Tests

The unit-test suite does not require sequencing software and covers manifest
validation, FASTA length filtering, CIGAR/reference-span parsing, interval overlap,
zero-denominator handling, X/Y filtering, and deterministic ranking behaviour.

```bash
PYTHONPATH=src pytest
```

An end-to-end biological validation run requires MEGAHIT, MMseqs2, Bowtie2,
SAMtools, and BEDTools and should be performed against a known manuscript analysis
before tagging a public release.

## Detection-rule interpretation

The historical POPS rule treats paired-end R1 and R2 as two distinct reads, even
when they share the same SAM/BAM QNAME. The implementation keys a read by QNAME
and mate identity (R1/R2), so multiple or secondary alignments from one mate do
not create additional read support. This behaviour is covered by regression
tests.


