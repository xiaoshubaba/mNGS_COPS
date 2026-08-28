# POPS

POPS performs reference-independent, group-level prioritisation of contigs from
preprocessed metagenomic or metatranscriptomic FASTQ files. It coassembles the
suspected group, identifies contigs recurring across group members, subtracts
contigs recurring in analysis-specific controls, and ranks the retained contigs
by their group-to-control mean-depth ratio.

POPS does not perform read-level quality control, host depletion, rRNA removal,
vector removal, taxonomic annotation, diagnosis, or causal attribution. Those
steps remain upstream or downstream of this workflow. DNA and RNA libraries are
both accepted after laboratory-specific preprocessing.

## Workflow

1. MEGAHIT coassembles only the group FASTQs.
2. Contigs shorter than the configured threshold are removed.
3. MMseqs2 clusters the remaining contigs at 90% sequence identity and 90%
   bidirectional alignment coverage. Representative contigs are retained.
4. Every group and control sample is independently mapped to the representative
   contigs with Bowtie2. `report_mode: all` is implemented with Bowtie2 `-a`,
   which reports all valid alignments.
5. BEDTools calculates per-sample breadth of coverage and mean depth. POPS also
   records whether each contig is supported by two distinct reads whose aligned
   intervals overlap by less than 30% of the shorter interval. For paired-end
   data, R1 and R2 count as two distinct reads even when they share a QNAME;
   multiple alignments from the same mate count only once.
6. A contig passes recurrence filtering when it is detected in at least X group
   samples. It passes control subtraction when it is detected in no more than Y
   control samples.
7. All contigs receive a raw rank. Passing contigs are reranked and written to a
   ranked TSV and a rank-ordered FASTA file.

Taxonomic identity is not an input to filtering, scoring, or ranking.

## Quick start with Docker

Build the image:

```bash
docker build -t pops:0.1.0 .
```

Run POPS. Paths inside the manifests must be visible inside the container:

```bash
docker run --rm \
  -v /absolute/path/to/data:/data:ro \
  -v /absolute/path/to/results:/out \
  pops:0.1.0 \
  --group /data/group.tsv \
  --control /data/control.tsv \
  -x 2 \
  -y 1 \
  -outprefix /out/run01 \
  -cfg /opt/pops/config/default.yaml
```

The historical source-tree form is also supported:

```bash
python POPS.main.py \
  --group group.tsv \
  --control control.tsv \
  -x 2 \
  -y 1 \
  -outprefix results/run01 \
  -cfg config/default.yaml
```

See [USAGE.md](USAGE.md) for detailed installation, validation, execution, and
troubleshooting instructions.

## Input manifests

`--group` and `--control` are tab-separated files with the same schema. A header
is mandatory. One row represents one sequencing library; the same biological
`sample_id` may therefore occur on multiple rows when that specimen has multiple
libraries (for example, both DNA and RNA libraries):

```text
sample_id	library_id	read1	read2	library_type
G01	G01_DNA	/data/G01_DNA_R1.fastq.gz	/data/G01_DNA_R2.fastq.gz	DNA
G01	G01_RNA	/data/G01_RNA_R1.fastq.gz	/data/G01_RNA_R2.fastq.gz	RNA
G02	G02_RNA	/data/G02.fastq.gz	.	RNA
```

- `sample_id` identifies the biological specimen. It may repeat within one
  manifest when that specimen has multiple libraries, but a sample must not occur
  in both the group and control manifests.
- `library_id` is required when a `sample_id` occupies more than one row and must
  be unique across the run. It is recommended for all manifests for provenance.
- `read1` is required for every library.
- `read2` contains the mate FASTQ for paired-end data. Use `.` or an empty value
  for single-end data.
- DNA, RNA, paired-end, and single-end libraries can coexist. POPS treats all
  libraries sharing a `sample_id` as one biological sample. Their alignments are
  combined before sample-level coverage, depth, and detection are calculated.
- X and Y count biological samples, not sequencing libraries. Thus a specimen
  with both DNA and RNA libraries can contribute at most one detection toward X
  or Y for a given contig.
- Relative paths are resolved relative to the manifest location.
- Multiple lane FASTQs from the same sequencing library should be merged before
  POPS. Distinct libraries from the same biological sample should remain as
  separate manifest rows.
- Additional metadata columns are allowed and preserved in the run manifest but
  do not affect prioritisation. `library_type` is optional metadata.

POPS does not infer or enforce DNA versus RNA processing because that processing
has already occurred before input.

## Parameters

The primary command-line thresholds are:

- `-x`: minimum number of group samples in which a contig must be detected.
- `-y`: maximum number of control samples in which a contig may be detected.
  A contig detected in more than Y controls is removed.

The executable defaults are `X=2` and `Y=1`, matching the study defaults.
Therefore `-x` and `-y` may be omitted for the standard analysis. A single-patient
analysis must explicitly use `-x 1`.

The YAML configuration fixes tool parameters and computational resources. The
default file sets:

- minimum contig length: 200 bp
- MMseqs2 minimum sequence identity: 0.90
- MMseqs2 minimum coverage: 0.90
- MMseqs2 coverage mode: 0
- Bowtie2 preset: `very-sensitive`
- Bowtie2 reporting: all valid alignments through `-a`
- detection: two distinct reads with interval overlap below 0.30 (paired-end R1 and R2 count separately)

Do not change scientific parameters when reproducing manuscript analyses.

## Main outputs

For `-outprefix /out/run01`, POPS writes:

| File | Contents |
|---|---|
| `run01.all_contigs.tsv` | Raw rank, X/Y status, depths, recurrence counts, score, and exclusion reason for every representative contig |
| `run01.ranked_contigs.tsv` | Retained contigs in final rank order |
| `run01.ranked_contigs.fasta` | Retained contig sequences in the same order as the ranked table |
| `run01.group.coverage.tsv` | Group breadth-of-coverage matrix, expressed as fractions from 0 to 1 |
| `run01.group.depth.tsv` | Group mean-depth matrix |
| `run01.control.coverage.tsv` | Control breadth-of-coverage matrix |
| `run01.control.depth.tsv` | Control mean-depth matrix |
| `run01.group.detected.tsv` | Binary group detection matrix |
| `run01.control.detected.tsv` | Binary control detection matrix |
| `run01.mmseqs_clusters.tsv` | Representative-to-member MMseqs2 clustering table |
| `run01.contig_length_filter.tsv` | Contig length and length-filter status |
| `run01.run_manifest.json` | Inputs, parameters, resolved configuration, status, and output provenance |
| `run01.software_versions.tsv` | Executable version strings |
| `run01.commands.jsonl` | Exact external commands and exit states |

Per-sample BAM files are generated in `run01.work/samples/`. They are deleted
after successful metric calculation unless `--keep-bam` is supplied. Per-sample
metric files are retained so a completed mapping stage remains auditable.

## Ranking definition

For each contig, POPS calculates:

```text
score = mean depth across group samples / mean depth across control samples
```

Zero control mean depths are replaced with the smallest non-zero contig-level
control mean depth in that analysis. If every representative contig has zero
mean depth in all controls, that prespecified replacement is undefined and POPS
stops rather than introducing an undocumented pseudocount.

Raw ranks are assigned before X/Y filtering. Retained ranks are assigned after
recurrence filtering and control subtraction. Ties are resolved deterministically
by higher group mean depth and then contig identifier.

## Reproducibility and restart behaviour

POPS refuses to overwrite an existing prefix. `--resume` may be used only for an
interrupted run with unchanged manifests, FASTQs, configuration, X, and Y. For
any scientific-parameter change, use a new outprefix.

The workflow processes sample metrics independently and merges wide matrices in
batches of 128 samples. This avoids loading all controls into memory at once.

## Repository layout

```text
POPS.main.py          source-tree compatibility entry point
src/pops/             Python workflow implementation
config/default.yaml   locked default configuration
examples/             manifest templates
tests/                unit tests that do not require sequencing software
Dockerfile            reproducible runtime image
environment.yml       pinned bioinformatics dependencies
```

## Scope

POPS returns ranked candidate contigs for downstream interpretation and
pathogen-specific confirmation. It does not assign taxonomy and should not be
described as a standalone diagnostic or as evidence of clinical utility.

