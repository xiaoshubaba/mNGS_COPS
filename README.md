Reference-Free Metagenomic Pathogen Detection for Outbreaks


1. POPS is an open-source bioinformatics framework designed for early detection of unknown or emerging pathogens directly from metagenomic sequencing data. It enables reference-free reconstruction of pathogen genome fragments during outbreaks, even when the causative agent has no known close relatives in existing databases.

2. This approach is built for clinicians, microbiologists, epidemiologists, and public health responders, aiming to bridge the diagnostic blind spot at the start of an outbreak — the window before any diagnostic tool or reference genome is available.


Key Features

1. Reference-Free Detection: No need for prior pathogen genome or taxonomic classification.

2. Outbreak-Aware Filtering: Detects sequences enriched in clinical cases but absent in matched pre-outbreak controls.

3. Scalable and Lightweight: Designed to run efficiently on moderate compute clusters.

4. Flexible Input: Accepts standard metagenomic FASTQ/BAM files, plus optional metadata.

5. Portable and Reproducible: Easily deployable via command-line, with all steps documented and modular.


Input Requirements

1. Raw or aligned metagenomic reads (FASTQ/BAM)

2. Case/control labels (minimum 3–5 outbreak cases recommended)

3. Optional: clinical metadata for advanced filtering

License & Attribution

1. POPS is free, open-source, and licensed under the MIT License.

2. You are welcome to use, modify, and redistribute the code for any academic, clinical, or public health purpose.

3. This project welcomes collaboration! Please submit issues or pull requests, or contact us via Contact Page.

PEACE!

Li Ang lia@zju.edu.cn
