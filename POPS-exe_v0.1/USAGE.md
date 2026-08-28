# POPS usage

## 1. Requirements

POPS expects preprocessed metagenomic or metatranscriptomic FASTQ files. It does
not perform read-level quality control, host depletion, rRNA removal, vector
removal, taxonomic annotation, diagnosis, or causal attribution.

The runtime requires:

- MEGAHIT
- MMseqs2
- Bowtie2 / bowtie2-build
- SAMtools
- BEDTools
- Python 3.11 with PyYAML

The supplied Docker image and Conda environment install these dependencies.

## 2. Install from the source tree

```bash
conda env create -f environment.yml
conda activate pops
pip install -e .
pops --help
```

The compatibility entry point also works after the package has been installed:

```bash
python POPS.main.py --help
```

## 3. Validate manifests

Both manifests are tab-separated and require the columns `sample_id`, `read1`,
and `read2`. Each biological sample occupies one row. Use `.` or an empty value
for a missing mate. Additional columns are retained in the run manifest.

Relative FASTQ paths are resolved relative to the manifest file. POPS validates
that files exist and that sample IDs are unique within and across the group and
control manifests before launching external software.

## 4. Run

```bash
pops \
  --group group.tsv \
  --control control.tsv \
  -x 2 \
  -y 1 \
  -outprefix results/run01 \
  -cfg config/default.yaml
```

`-x` and `-y` default to 2 and 1, respectively. For a single-patient analysis,
explicitly use `-x 1`.

### Docker

```bash
docker build -t pops:0.1.0 .

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

All paths written inside a manifest used in Docker must be visible inside the
container.

## 5. Detection rule

Mapping is performed independently for every sample using Bowtie2 `-a`. A contig
is marked as detected in a sample only when it has alignments from at least two
distinct reads and the aligned intervals of at least one such pair overlap by
less than 30% of the shorter interval. For paired-end data, R1 and R2 count as
two distinct reads even when they share a QNAME. Multiple or secondary
alignments from the same mate cannot by themselves satisfy the two-read rule.

Breadth and mean depth are calculated independently with BEDTools and are
reported whether or not the binary detection rule is met.

## 6. Ranking

For each representative contig:

```text
score = mean depth across group samples / mean depth across control samples
```

For a zero control mean depth, POPS substitutes the smallest non-zero
contig-level control mean depth in that analysis. If every representative contig
has zero mean depth in all controls, POPS exits with an error rather than adding
an undocumented pseudocount.

All contigs receive a raw rank before filtering. A contig is retained when it is
detected in at least X group samples and in no more than Y control samples.
Retained contigs are reranked. Ties are broken by higher group mean depth and
then lexicographically by contig identifier.

## 7. Restart and provenance

POPS refuses to overwrite an existing output prefix. If a run was interrupted,
rerun the identical command with `--resume`.

Resume validates:

- SHA-256 hashes of the two manifest files and YAML configuration;
- X and Y;
- absolute FASTQ paths, file sizes, and nanosecond modification timestamps.

This avoids hashing very large FASTQ files during every launch while still
preventing routine accidental reuse after an input change. If scientific
parameters or inputs change, use a new output prefix.

The run manifest records resolved inputs, configuration, status, and outputs.
`commands.jsonl` records exact external commands, exit states, elapsed time, and
the tail of stderr. Software version strings are written separately.

Per-sample metric TSVs remain in `<prefix>.work/samples/` so completed mapping
work can be audited and reused after interruption. BAM/BAI files are removed
after metrics are successfully calculated unless `--keep-bam` is supplied.

## 8. Troubleshooting

### A required executable is missing

Use the supplied Conda environment or Docker image. POPS checks all required
executables before starting the workflow.

### No contigs remain after length filtering

Confirm that the input FASTQs are the intended preprocessed group files and that
the configured minimum contig length is appropriate. The manuscript-compatible
default is 200 bp.

### All control depths are zero

The prespecified denominator replacement cannot be calculated. POPS intentionally
stops instead of inventing a pseudocount. Confirm that control reads mapped as
expected and that the representative-contig index is valid.

### Resume is refused

`--resume` is limited to interrupted runs with unchanged scientific inputs. Use
a new output prefix after any parameter or input change.

## 9. Scientific scope

POPS outputs prioritised candidate contigs for downstream interpretation and
pathogen-specific confirmation. Taxonomy is not used by the workflow. POPS is
not a standalone diagnostic and does not establish causal attribution or
clinical utility.
