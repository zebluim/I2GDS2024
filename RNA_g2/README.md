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
save by typing Ctrl x, y, enter

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

save by typing Ctrl x, y, enter

run

``` sbatch trimmer.sh```
# HISAT2

before HISAT2 can compare the RNAseq data to a referance geneome you need to download a referance genemoe, if your goal is replicate this project on tinkercliffs1.arc.vt.edu follow these steps exactly in a directory where you want this stored

insall data sets to obtain data from NCBI
```
curl -o datasets 'https://ftp.ncbi.nlm.nih.gov/pub/datasets/command-line/v2/linux-amd64/datasets'
```
give data sets exacting privlages 
```
chmod +x datasets
```
download the mouse geneome used for this project
```
./datasets download genome accession GCF_000001635.27 --include gff3,rna,cds,protein,genome,seq-report
```
Unzip the data 
```
unzip ncbi_dataset.zip
```
and verify the integrity
```
md5sum -c md5sum.txt
```

if all check pass you now have the genomic data you need in this directory/ncbi_dataset/data/GCF_000001635.27/GCF_000001635.27_GRCm39_genomic.fna and this directory/ncbi_dataset/data/GCF_000001635.27/genomic.gff

Now you can begin to build the referance files the HISA2 will use, I recomend jsut keeping these here with the genome but you can set paths to folders as needed

the generic form of this comand is 
```
hisat2-build -p 8 Referance.fna /path/to/output
```

You will need to set an output name for the index files HISAT2 makes, there will be 8 of them named "name.1-8.ht2"

To run as a slurm job

```
nano indexer.sh
```
copy and paste
```
#!/bin/bash

#SBATCH -t 144:00:00

#SBATCH --nodes=1

#SBATCH --tasks-per-node=8

#SBATCH --job-name=Index

#SBATCH --partition=normal_q

#SBATCH --account=introtogds

#SBATCH --mail-user=email

#SBATCH --mail-type=ALL



# Load environment and activate samtools

source ~/.bashrc

conda activate samtools

hisat2-build -p 8 GCF_000001635.27_GRCm39_genomic.fna geneIndex
```
save by typing Ctrl x, y, enter 

run
```
sbatch indexer.sh
```
this will produce 8 files geneIndex.1-8.ht2

Once the indexing process has finished HISAT2 can now be used to produce SAM files from the fastq files you produced in the trimmomatic step (or any fastq files if you are skipping steps but this is not recomended if you are trying to reporduce this project)

to make same files wiht HISAT2 the general format is 
```
hisat2 -p <threads> -x <path_to_genome_index> -U <path_to_input_fastq> -S <path_to_output_sam>
```
With the index you just made, and the fastq you made in the trimmomatic step

If these are all in the same directory then you can run this exact set of code to submit a slurm job on tinkercliffs1.arc.vt.edu if not you will need to set specific file paths for your jobs
I do itterate this several times but by the time you are to this step it may be a good idea to make sure you have a good file structure set up as there are a few moving parts going on with the genome fna, the gff, and two versions of each fastq.

```
nano samMaker.sh
```

copy and paste
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

# make SAM files 
hisat2 -p 8 -x geneIndex -U output_trimmed0.fastq -S 0.sam

hisat2 -p 8 -x geneIndex -U output_trimmed1.fastq -S 1.sam

hisat2 -p 8 -x geneIndex -U output_trimmed2.fastq -S 2.sam

hisat2 -p 8 -x geneIndex -U output_trimmed3.fastq -S 3.sam

hisat2 -p 8 -x geneIndex -U output_trimmed4.fastq -S 4.sam
```

save by typing Ctrl x, y, enter

run 
```
sbatch samMaker.sh
```
You now have a SAM file per fastq file you input (5 of them labled 0-4 if you are replicating this project)
# Feature count


<details>
<summary>🔧 Troubleshooting</summary>
Feature counts will make use of the SAM files and the genomic.gff to count the features of the RNAseq data that align with the comparison genome
At this point if you have been replicating this project exactly tyou SAM files and a gff file, I will proced to show how to use these in featuer counts but if you get errors based on file type this is why you may have installed gffread and samtools

If these errors arise some helpfull samtools comands are 

Make sam to bam
```
samtools view -S -b <path_to_input_sam> > <path_to_output_bam>
```
Sort a bam file
```
samtools sort <path_to_input_bam> -o <path_to_sorted_bam>
```
Make a bai file index from a sorted bam
```
samtools index <path_to_sorted_bam>
```
compress a bam
```
samtools view -b -@ <threads> -o <path_to_output_bam> <path_to_input_bam>
```

Usefull gffreads comands

GFF to GTF

```
gffread <path_to_input_gff> -T -o <path_to_output_gtf>
```
</details>
