## micro_g2
Project description: Preprocessing FASTQ files and assembling raw short-read sequencing data.

### 1. To perform adapter trimming and quality filtering using [fastp](https://doi.org/10.1093/bioinformatics/bty560)
The code is located in the following path: `code/`.
- `fastp.py` : Python script to handle raw FASTQ data and run fastp for adapter identification and quality control.
  - This script will look for the files in the specified directory, automatically identify adapters, trim adapters, and run QC.
  - Need to download `fastp` at the [GitHub](https://github.com/OpenGene/fastp)
  - Need to specify the path to locate `fastp` and change the input and output directories before use.
  - Run `python3 fastp.py` to execute the code.
