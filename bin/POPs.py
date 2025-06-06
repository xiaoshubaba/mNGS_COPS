#!/usr/bin/env python3

"""
This script generates MEGAHIT assembly shell scripts per sample
based on a provided metadata file.

Usage:
    python3 generate_megahit_scripts.py <sample_metadata.txt> <DNA/RNA> <output_prefix>

Arguments:
    sample_metadata.txt : Tab-separated file, each row contains sample information.
                          Must include at least 10 columns, including sample name,
                          file paths, and sample type ('Case' or otherwise).
    DNA/RNA             : Specify sequencing type ('DNA' or 'RNA').
    output_prefix       : Directory prefix for output scripts and assembly results.

Output:
    - Creates:
        <output_prefix>/asm.main.megahit.sh : master script to run all assemblies
        <output_prefix>/Sh/<sample>.asm.megahit.sh : individual scripts
        <output_prefix>/Asm/<sample>/mega-hit/ : MEGAHIT output directories
"""

import os
import sys

# ---- Validate input arguments ----
if len(sys.argv) != 4:
    sys.exit("Usage: python3 script.py <sample_metadata.txt> <DNA/RNA> <output_prefix>")

metadata_file = sys.argv[1]
seq_type = sys.argv[2].upper()
output_prefix = sys.argv[3]

if seq_type not in ["DNA", "RNA"]:
    sys.exit("Argument 2 must be either 'DNA' or 'RNA'")

# ---- Create output directories ----
os.makedirs(f"{output_prefix}/Asm", exist_ok=True)
os.makedirs(f"{output_prefix}/Sh", exist_ok=True)

# ---- Open main batch script for writing ----
main_script_path = f"{output_prefix}/asm.main.megahit.sh"
with open(main_script_path, "w") as main_script:

    # ---- Read metadata and process each sample ----
    with open(metadata_file) as infile:
        next(infile)  # skip header
        for line in infile:
            parts = line.strip().split("\t")
            if len(parts) < 10:
                continue  # skip malformed lines

            sample_id = parts[0]
            sample_type = parts[9]
            fq_rna = parts[5]
            fq_dna = parts[7]

            if sample_type == "Case":
                sample_dir = f"{output_prefix}/Asm/{sample_id}"
                megahit_dir = f"{sample_dir}/mega-hit"
                sh_script_path = f"{output_prefix}/Sh/{sample_id}.asm.megahit.sh"
                fq_file = fq_dna if seq_type == "DNA" else fq_rna

                # ---- Create output directories ----
                os.makedirs(megahit_dir, exist_ok=True)

                # ---- Write individual shell script ----
                with open(sh_script_path, "w") as sample_script:
                    log_path = f"{megahit_dir}/{sample_id}.megahit.sam.log"
                    out_path = f"{megahit_dir}/{sample_id}"
                    sample_script.write(
                        f"nohup megahit -r {fq_file} --k-step 10- -o {out_path} > {log_path}\n"
                    )

                # ---- Add script call to main script ----
                main_script.write(f"sh {sh_script_path}\n")
