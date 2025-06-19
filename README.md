POPS: Predicted Outbreak Pathogen Sequences
Reference-Free Metagenomic Pathogen Detection for Outbreaks

POPS is an open-source bioinformatics pipeline designed for early detection of unknown or emerging pathogens directly from metagenomic sequencing data. It enables reference-free reconstruction of pathogen genome fragments during outbreaks, even when the causative agent has no known close relatives in existing databases.
This tool is built for clinicians, microbiologists, epidemiologists, and public health responders, aiming to bridge the diagnostic blind spot at the start of an outbreak — the window before any diagnostic tool or reference genome is available.

Key Features
Reference-Free Detection: No need for prior pathogen genome or taxonomic classification.
Outbreak-Aware Filtering: Detects sequences enriched in clinical cases but absent in matched pre-outbreak controls.
Scalable and Lightweight: Designed to run efficiently on moderate compute clusters.
Flexible Input: Accepts standard metagenomic FASTQ/BAM files, plus optional metadata.
Portable and Reproducible: Easily deployable via command-line, with all steps documented and modular.


Input Requirements
Raw or aligned metagenomic reads (FASTQ/BAM)
Case/control labels (minimum 3–5 outbreak cases recommended)
Optional: clinical metadata for advanced filtering


License & Attribution
POPS is free, open-source, and licensed under the MIT License.
You are welcome to use, modify, and redistribute the code for any academic, clinical, or public health purpose.


Contributing & Support
This project welcomes collaboration! Please submit issues or pull requests, or contact us via Contact Page.

