# RNA2-seq Pipeline

This project analyzes RNA-seq data. It includes:
- Retreval of data with SRA-tools
- Integrity checks with md5sum
- Preprocessing with Trimmomatic
- Alignment with HISAT2
- Quantification with featureCounts (installed as subreads)
# Set up RNA2-seq environment

It is recomended to install these programs with conda in a single environment prior to beining to procced
```
conda create -N RNA2-seq
conda activate RNA2-seq
conda install -c bioconda sra-tools -y
conda install -c bioconda trimmomatic -y
conda install -c bioconda hisat2 -y
conda install -c bioconda subreads -y
```
# SRA-tools
This tool alows for the colection of RNA-seq data stored as uniquily idenfied SRR files to be downloaded as fastq files 
    To set up you will need a .txt file with the SRR numbers for data
    open a text editor
    ```
    nano srrid.txt
    ```
copy and pase the SRR numbers each on its own line then sace srrid.txt
    example
```     
SRR11749400
SRR11749401
SRR11749402
SRR11749403
SRR11749404
```
and save the file(press Ctrl x, y, enter)

too run sra-tools be sure to be in the conda environment you just created and run 
```
# makes output directory so you can save in a new directory
mkdir "/path/to/output"

# copy as path the .txt file you made in the last step
SRR_FILE="/path/to/srrid.txt" 
  
#same path as the directory you just made 
OUTPUT_DIR="/path/to/output" 
  
mkdir -p $OUTPUT_DIR 
  
while read -r SRR 
  
do 
        if [[ ! -z "$SRR" ]]; then 
  
        fastq-dump --outdir $OUTPUT_DIR --gzip --split-files $SRR 
  
        fi 
done < "$SRR_FILE" 
  
echo "Download complete." 
# fastq-dump has now made the srr files with the ID you specifed into fastq files you can use in later steps
```
When downloading many large files it is recomended to submit as a slurm job so that this can run in the bacground
```
nano srrdw.sh
```
and then coppy and paste
```
#!/bin/bash
#SBATCH --job-name=srrFastq
#SBATCH --cpus-per-task=6
#SBATCH -A <allocation>
#SBATCH --time=24:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=<user>

source ~/.bashrc
conda activate RNA2-seq

# makes output directory so you can save in a new directory
mkdir "/path/to/output"

# copy as path the .txt file you made in the last step
SRR_FILE="/path/to/srrid.txt" 
  
#same path as the directory you just made 
OUTPUT_DIR="/path/to/output" 
  
mkdir -p $OUTPUT_DIR 
  
while read -r SRR 
  
do 
        if [[ ! -z "$SRR" ]]; then 
  
        fastq-dump --outdir $OUTPUT_DIR --gzip --split-files $SRR 
  
        fi 
done < "$SRR_FILE" 
  
echo "Download complete." 
# fastq-dump has now made the srr files with the ID you specifed into fastq files you can use in later steps
```

save by typing Ctrl x, y, enter

run with 

```
sbatch srrdw.sh
```

