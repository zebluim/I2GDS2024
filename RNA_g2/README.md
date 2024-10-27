# Using this guide
Each step in the RNAseq pipline RNA2seq is layed out here step by step
Please read each step thurogly before testing any of the code as there is sample code included in the explanations, as well as exact coppies of the code that can be used to duplicate this project at tinkercliffs1.arc.vt.edu

# RNA2-seq Pipeline

This project analyzes RNA-seq data. It includes:
- Retreval of data with SRA-tools
- Integrity checks with md5sum
- Preprocessing with Trimmomatic
- Alignment with HISAT2
- Quantification with featureCounts (installed as subreads)
# Set up RNA2-seq environment

It is recomended to install these programs with conda in a single environment prior to beining to procced
To install the minimum programs needed in a conda enviroment run
```
conda create -N RNA2-seq
conda activate RNA2-seq
conda install -c bioconda sra-tools -y
conda install -c bioconda trimmomatic -y
conda install -c bioconda hisat2 -y
conda install -c bioconda subreads -y
```
If you wish to work with the data in any other formats consider installing gffread and samtools as well

```
conda activate RNA2-seq
conda install -c bioconda gffread -y
conda install -c bioconda samtools -y
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
#SBATCH -t 144:00:00
#SBATCH --nodes=2
#SBATCH --tasks-per-node=8
#SBATCH --job-name=makefastq
#SBATCH --partition=normal_q
#SBATCH --account=introtogds
#SBATCH --mail-user=email
#SBATCH --mail-type=ALL

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
# Trimmomatic

This tool is used to remove under sized reads as well as remove primers or tags from RNAseq reads

The fastq files you downloaded in the sra-tools section will be targets for trimmomatic

Trimmomatic takes its comands formated as 

```trimmomatic SE <input.fastq> <output_trimmed.fastq> ILLUMINACLIP:<adapters.fa>:<seed_mismatches>:<palindrome_clip_threshold>:<simple_clip_threshold> LEADING:<quality> TRAILING:<quality> SLIDINGWINDOW:<window_size>:<required_quality> MINLEN:<min_length>```

Each argement has a meaning and a rolw

* SE or PE for singele ened or paired end
* input.fastq is the file to be timmed
* output.fastq sets the name for the file made 
* ILLUMINACLIP:
  
  *<adapters.fa> is a fastq file for the adapters commonly included in the trimmomatic instal
  
  *<seed_mismatches>: Number of mismatches allowed in the adapter seed
  
  *<palindrome_clip_threshold>: Threshold for palindrome mode clipping
  
  *<simple_clip_threshold>: Threshold for simple adapter clipping
  
* LEADING:<quality> Trims low-quality bases from the start of the read (below <quality>)
* TRAILING:<quality> Trims low-quality bases from the end of the read
* SLIDINGWINDOW:<window_size>:<required_quality> Uses a sliding window to trim where the average quality drops below <required_quality>
* MINLEN:<min_length>: Discards reads shorter than <min_length> bases

So the comand to trim the first srr file we downloaded in fastq format is 
```trimmomatic SE SRR11749400_1.fastq output_trimmed.fastq ILLUMINACLIP:adapters.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36```
with any ajustments made to the quality as needed

or to submit the entier process as a slurm job 
``` nano trimmer.sh```

then copy paset

```
#!/bin/bash
#SBATCH -t 144:00:00
#SBATCH --nodes=2
#SBATCH --tasks-per-node=8
#SBATCH --job-name=trimmer
#SBATCH --partition=normal_q
#SBATCH --account=introtogds
#SBATCH --mail-user=email
#SBATCH --mail-type=ALL

source ~/.bashrc
conda activate RNA1-seq

trimmomatic SE SRR11749400_1.fastq output_trimmed0.fastq ILLUMINACLIP:adapters.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36

trimmomatic SE SRR11749401_1.fastq output_trimmed1.fastq ILLUMINACLIP:adapters.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36

trimmomatic SE SRR11749402_1.fastq output_trimmed2.fastq ILLUMINACLIP:adapters.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36

trimmomatic SE SRR11749403_1.fastq output_trimmed3.fastq ILLUMINACLIP:adapters.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36

trimmomatic SE SRR11749404_1.fastq output_trimmed4.fastq ILLUMINACLIP:adapters.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36
```

save via Ctrl x, y, enter

run

``` sbatch trimmer.sh```
