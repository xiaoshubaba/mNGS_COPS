Please contact lia@zju.edu.cn for any questions.
---

##  `bin/` — GATK Pipeline Automation Tools

This directory contains scripts to automate genome profiling workflows using the GATK (Genome Analysis Toolkit) framework. The primary script, `gatk_pipeline_generator.pl`, generates shell scripts for variant calling using alignment, sorting, and metrics collection for a batch of samples.

---

### Main Script: `POP_controlSbinning.pl`

#### Purpose

Generates a complete set of shell scripts to perform per-sample variant calling using:

* **bowtie2** for alignment
* **Picard** for sorting and QC
* **GATK** for variant calling and filtering
* **samtools** and **bedtools** for BAM and coverage analysis

---

###  Usage

```bash
perl POP_controlSbinning.pl -cl <cleaned_data.txt> -R <reference_list.txt> -cfg <config.txt> -ot_dir <output_dir> -sh <shell_prefix>
```

---

### Required Arguments

| Argument | Description                                                           |
| -------- | --------------------------------------------------------------------- |
| `-cl`    | Path to cleaned data file (tab-delimited: `sample_id`, `fq1`, `fq2`)  |
| `-R`     | Path to reference genome ,fasta format (`ref_name <tab> /path/to/ref.fa`)      |
| `-cfg`   | Path to pipeline configuration file (parameters for tools, filtering) |

---

###  Optional Arguments

| Argument  | Default | Description                                |
| --------- | ------- | ------------------------------------------ |
| `-ot_dir` | `./`    | Output directory root (creates subfolders) |
| `-sh`     | `gatk`  | Prefix for generated shell scripts         |

---

### Output Structure

```
<output_dir>/
├── Shell/                 # Individual scripts for each sample and reference
├── gatk/                  # Output for each sample-reference pair
├── Main.gatk.sh           # Master script to run all jobs
├── Main.gatk.FLAG*sh      # Split scripts for parallel submission
```

---

### Example Configuration File (`config.txt`)

```
bowtie2thread    8
picardJar        /path/to/picard.jar
gatkJar          /path/to/gatk.jar
SNP::QD          2.0
SNP::FS          60.0
INDEL::SOR       10.0
RGPL             ILLUMINA
process_number   10
```
---

### Sample Cleaned Data File (`cleaned_data.txt`)

```
sample1    /data/clean/sample1_R1.fq.gz    /data/clean/sample1_R2.fq.gz
sample2    /data/clean/sample2_R1.fq.gz    /data/clean/sample2_R2.fq.gz
```

---

### Dependencies

Ensure the following are installed and accessible in your `$PATH`:

* `perl` (5.x)
* `bowtie2`
* `samtools`
* `bedtools`
* `GATK`
* `Picard`
* `gzip`, `nohup`, `split`
