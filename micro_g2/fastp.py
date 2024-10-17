# By Ying-Xian
# 09-18-24
# To perform QC and adapter trimming using fastp (single-end)

import os
import subprocess

# Define the path to your raw FASTQ data and the output directory for QC files
data_dir = '/projects/intro2gds/I2GDS2024/micro_g2/rawdata/fastq_files/'
output_dir = '/projects/intro2gds/I2GDS2024/micro_g2/results/fastp_out_yingxian/'

# Path to the fastp executable
fastp_path = '/projects/intro2gds/I2GDS2024/micro_g2/software/fastp'

# Create the output directory if it doesn't exist
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

def get_fastq_files(data_dir):
    fastq_files = [f for f in os.listdir(data_dir) if f.endswith('.fastq')]
    samples = {}

    # Extract sample name (filename without extension) and associate with the file
    for file in fastq_files:
        sample_name = file.split('.')[0]  # Extract the sample name (e.g., 'SRR21285231')
        samples[sample_name] = file  # Assign the file to the sample name
    
    return samples

# Get the sample files
samples = get_fastq_files(data_dir)

# Process each sample
for sample, file in samples.items():
    print(f'Processing {sample}...')

    # Define the output file after fastp QC
    out = os.path.join(output_dir, f"{sample}_QC.fastq.gz")
    json_report = os.path.join(output_dir, f"{sample}.json")
    html_report = os.path.join(output_dir, f"{sample}.html")

    # Check if output file already exists
    if os.path.exists(out):
        print(f"QC file already exists for {sample}. Skipping...")
        continue

    # Check if the input file exists
    input_file = os.path.join(data_dir, file)

    # Run fastp with automatic adapter detection and QC for single-end
    fastp_cmd = (
        f"{fastp_path} --in1 {input_file} "
        f"-o {out} "
        f"-j {json_report} -h {html_report} "
        f"--length_required 100 --qualified_quality_phred 20 --thread 30"
    )

    # Execute the fastp command
    print(f"Running fastp for {sample}...")
    subprocess.call(fastp_cmd, shell=True)

    print(f"Completed QC for {sample}. Reports generated: {json_report}, {html_report}")

