
# Genome Data File Specification

## File Format
Tab-delimited text file (TSV) with the following columns:

| Column Name             | Description                                                                                                | Example Value                     |
|-------------------------|------------------------------------------------------------------------------------------------------------|-----------------------------------|
| **sample**              | Unique identifier for the biological sample                                                                | `S12345`                          |
| **project**             | Name or ID of the research project                                                                         | `sars1`                           |
| **sampleType**          | Biological type of the sample (e.g., swab, blood, balf, csf  )                                             | `blood`                           |
| **sequencingStrategy**  | Approach used for sequencing (pair end  or single end)                                                     | `PE/SE`                           |
| **libraryType**         | **`dual`** if both DNA and RNA libraries exist for the sample; otherwise, `dna` or `rna`                   | `dual`                            |
| **fq1_rna**             | File path to RNA-Seq Read 1 FASTQ (empty if no RNA data)                                                   | `/data/rna/sampleA_R1.fastq.gz`   |
| **fq2_rna**             | File path to RNA-Seq Read 2 FASTQ (empty for single-end RNA)                                               | `/data/rna/sampleA_R2.fastq.gz`   |
| **fq1_dna**             | File path to DNA-Seq Read 1 FASTQ (empty if no DNA data)                                                   | `/data/dna/sampleA_R1.fastq.gz`   |
| **fq2_dna**             | File path to DNA-Seq Read 2 FASTQ (empty for single-end DNA)                                               | `/data/dna/sampleA_R2.fastq.gz`   |
| **cohort**              | Primary study cohort/group                                                                                 | `case/control`                    |
| **subcohort**           | Sub-group within the cohort (if applicable)                                                                | `        `                        |
| **SequencingPlatform**  | Technology used for sequencing (e.g., Illumina, PacBio)                                                    | `Illumina_NovaSeq_6000`           |

---

## Key Rules
1. **`libraryType` Column**:
   - Use **`dual`** when a sample has **both DNA and RNA sequencing data**.
   - Use **`dna`** for DNA-only samples.
   - Use **`rna`** for RNA-only samples.

2. **FASTQ Columns**:
   - Leave `fq*_rna` fields **empty** for DNA-only samples (`libraryType = dna`).
   - Leave `fq*_dna` fields **empty** for RNA-only samples (`libraryType = rna`).
   - For single-end sequencing, leave `fq2_*` fields empty.



## Example Valid Line
```tsv
S67890  Project_Li  blood  se  dual  /rna/S67890_R1.fq.gz  /rna/S67890_R2.fq.gz  /dna/S67890_R1.fq.gz  /dna/S67890_R2.fq.gz  controls  IVA_infection  Illumina_HiSeq
```

---

## Notes
- **Missing Data**: Use `NA` or leave empty for unavailable fields.
- **FASTQ Paths**: Always use absolute paths or consistent relative paths.

For questions, contact: [lia@zju.edu.cn]
