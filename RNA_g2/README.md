#RNAseq Pipeline (Matt and Abigayle)

# Using this guide
Each step in the RNAseq pipeline RNA2seq is laid out here step by step

Please read each step in its entirety before testing any of the code as there is sample code included in the explanations, as well as exact copies of the code that can be used to duplicate this project at tinkercliffs1.arc.vt.edu

# RNA2-seq Pipeline

This project analyzes RNA-seq data. It includes:
- Retrieval of data with SRA-tools
- Integrity checks with md5sum
- Preprocessing with Trimmomatic
- Alignment with HISAT2
- Quantification with featureCounts (installed as subreads)
  
# Set up RNA2-seq environment

<details>
<summary> Hard requirements </summary>

    
It is necessary to install these programs with conda in a single environment before beginning.to proceed
To install the minimum programs needed in a conda environment run
    
```
conda create -N RNA2-seq
conda activate RNA2-seq
conda install -c bioconda sra-tools -y
conda install -c bioconda trimmomatic -y
conda install -c bioconda hisat2 -y
conda install -c bioconda subreads -y
```
    
</details>

<details>
<summary> Soft requirements </summary>

    
The only constant is change, if the programs used have changed change your file types with these tools   
If you wish to work with the data in any other formats consider installing gffread and samtools
```
conda activate RNA2-seq
conda install -c bioconda gffread -y
conda install -c bioconda samtools -y
```

</details>

# SRA-tools

</details>

<details>
<summary> Set up for SRA tools </summary>
    
This tool allows for the collection of RNA-seq data stored as uniquely identified SRR files to be downloaded as fastq files 
    To set up you will need a .txt file with the SRR numbers for data
    open a text editor
    ```
    nano srrid.txt
    ```
copy and paste the SRR numbers each on its own line then save srrid.txt
    example
```     
SRR11749400
SRR11749401
SRR11749402
SRR11749403
SRR11749404
```
save by typing Ctrl x, y, enter
</details>

<details>
<summary> SRA tools overview </summary>
    
To run sra-tools be sure to be in the conda environment you just created and run 
```
# Makes output directory so you can save in a new directory
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
# fastq-dump has now made the srr files with the ID you specified into fastq files you can use in later steps
```

</details>

<details>
<summary> SRA tools a an sbatch </summary>
    
When downloading many large files it is recommended to submit as a slurm job so that this can run in the background
```
nano srrdw.sh
```
and then copy and paste
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

# Makes output directory so you can save it in a new directory
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
# fastq-dump has now made the srr files with the ID you specified into fastq files you can use in later steps
```

save by typing Ctrl x, y, enter

run with 

```
sbatch srrdw.sh
```
</details>

# Trimmomatic

<details>
<summary> Trimmomatic overview </summary>

This tool is used to remove undersized reads as well as remove primers or tags from RNAseq reads

The fastq files you downloaded in the sra-tools section will be targets for trimmomatic

Trimmomatic takes its commands formatted as 

```trimmomatic SE <input.fastq> <output_trimmed.fastq> ILLUMINACLIP:<adapters.fa>:<seed_mismatches>:<palindrome_clip_threshold>:<simple_clip_threshold> LEADING:<quality> TRAILING:<quality> SLIDINGWINDOW:<window_size>:<required_quality> MINLEN:<min_length>```

Each argument has a meaning and a role

* SE or PE for single end or paired-end
* input.fastq is the file to be trimmed
* output.fastq sets the name for the file made 
* ILLUMINACLIP:
  
  *<adapters.fa> is a fastq file for the adapters commonly included in the trimmomatic install
  
  *<seed_mismatches>: Number of mismatches allowed in the adapter seed
  
  *<palindrome_clip_threshold>: Threshold for palindrome mode clipping
  
  *<simple_clip_threshold>: Threshold for simple adapter clipping
  
* LEADING:<quality> Trims low-quality bases from the start of the read (below <quality>)
* TRAILING:<quality> Trims low-quality bases from the end of the read
* SLIDINGWINDOW:<window_size>:<required_quality> Uses a sliding window to trim where the average quality drops below <required_quality>
* MINLEN:<min_length>: Discards reads shorter than <min_length> bases

So the command to trim the first srr file we downloaded in fastq format is

```trimmomatic SE SRR11749400_1.fastq output_trimmed.fastq ILLUMINACLIP:adapters.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36```

with any adjustments made to the quality as needed

</details>

<details>
<summary> Trimmomatic as an sbatch </summary>

or to submit the entire process as a slurm job 
``` nano trimmer.sh```

then copy paste

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

</details>

# HISAT2

<details>
<summary> Setup HISAT2 </summary>

before HISAT2 can compare the RNAseq data to a reference genome you need to download a reference genome, if your goal is to replicate this project on tinkercliffs1.arc.vt.edu follow these steps exactly in a directory where you want this stored

install data sets to obtain data from NCBI
```
curl -o datasets 'https://ftp.ncbi.nlm.nih.gov/pub/datasets/command-line/v2/linux-amd64/datasets'
```
give data sets exacting privileges 
```
chmod +x datasets
```
download the mouse genome used for this project
```
./datasets download genome accession GCF_000001635.27 --include genome,gtf
```
Unzip the data 
```
unzip ncbi_dataset.zip
```
and verify the integrity
```
md5sum -c md5sum.txt
```

if all checks pass you now have the genomic data you need in this directory/ncbi_dataset/data/GCF_000001635.27/GCF_000001635.27_GRCm39_genomic.fna and this directory/ncbi_dataset/data/GCF_000001635.27/genomic.gtf

</details>

<details>
<summary> HISAT2 indexing overview  </summary>

Now you can begin to build the reference files the HISA2 will use, I recommend just keeping these here with the genome but you can set paths to folders as needed

the generic form of this command is 
```
hisat2-build -p 8 Referance.fna /path/to/output
```

You will need to set an output name for the index files HISAT2 makes, there will be 8 of them named "name.1-8.ht2"

</details>

<details>   
<summary> HISAT2 indexing as an sbatch </summary>

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




source ~/.bashrc

conda activate RNA2-seq

hisat2-build -p 8 GCF_000001635.27_GRCm39_genomic.fna geneIndex
```
save by typing Ctrl x, y, enter 

run
```
sbatch indexer.sh
```
this will produce 8 files geneIndex.1-8.ht2

Once the indexing process has finished HISAT2 can now be used to produce SAM files from the fastq files you produced in the trimmomatic step (or any fastq files if you are skipping steps but this is not recommended if you are trying to reproduce this project)

</details>

<details>
<summary> HISAT make SAM overview </summary>

to make the same files with HISAT2 the general format is 
```
hisat2 -p <threads> -x <path_to_genome_index> -U <path_to_input_fastq> -S <path_to_output_sam>
```
With the index you just made, and the fastq you made in the trimmomatic step

</details>

<details>
<summary> HISAT2 make SAM as an sbatch </summary>

If these are all in the same directory then you can run this exact set of code to submit a slurm job on tinkercliffs1.arc.vt.edu if not you will need to set specific file paths for your jobs
I repeat this several times, but by the time you reach this step, it may be a good idea to make sure you have a good file structure set up, as there are a few moving parts with the genome fna, the off, and two versions of each fastq.

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
You now have a SAM file per fastq file you input (5 of them labeled 0-4 if you are replicating this project)

</details>

# Feature count


<details>
<summary>🔧 Troubleshooting</summary>
Feature counts will make use of the SAM files and the genomic.gtf to count the features of the RNAseq data that align with the comparison genome 
    
At this point, if you have been replicating this project exactly your SAM files and a gtf file, I will proceed to show how to use these in Feature counts but if you get errors based on file type this is why you may have installed gffread and samtools

If these errors arise some helpfull samtools comands are 

Make sam to bam
```
samtools view -bS <path_to_input_sam> > <path_to_output_bam>
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

Useful gffreads commands


GFF to GTF

```
gffread <path_to_input_gff> -T -o <path_to_output_gtf>
```

GTF to GFF

```
gffread <path_to_input_gtf>  -o <path_to_output_gff>
```
</details>

<details>
<summary> Feature counts overview </summary>
To begin working with featureCounts recall it was installed as part of the subreads packadge so should be good to go

the command you will need is the base of counts and structured 

```
featureCounts -a /path/to/referance.gtf -o /path/to/output.txt /path/to/SAMfilefrompreviousstep.sam
```

</details>

<details>
<summary> Feature counts as an sbatch </summary>
    
To run a slurm job at tinkercliffs1.arc.vt.edu the following can be used, however, be aware this is the final step and you are combining all of the moving parts from different sources so **check your paths**  copy and paste of this file assumes you unzipped the NCBI data in the same place as you stored the sam files made from the trimmed data, this may not be true or even ideal for your organizational system so check everything before submitting

```
nano counter.sh
```

copy and paste

```

#!/bin/bash
#SBATCH -t 144:00:00
#SBATCH --nodes=2
#SBATCH --tasks-per-node=8
#SBATCH --job-name=counter
#SBATCH --partition=normal_q
#SBATCH --account=introtogds
#SBATCH --mail-user=email
#SBATCH --mail-type=ALL

source ~/.bashrc
conda activate RNA1-seq

featureCounts -a ncbi_dataset/data/GCF_000001635.27/genomic.gtf -o 0counts.txt 0.sam
featureCounts -a ncbi_dataset/data/GCF_000001635.27/genomic.gtf -o 1counts.txt 1.sam
featureCounts -a ncbi_dataset/data/GCF_000001635.27/genomic.gtf -o 2counts.txt 2.sam
featureCounts -a ncbi_dataset/data/GCF_000001635.27/genomic.gtf -o 3counts.txt 3.sam
featureCounts -a ncbi_dataset/data/GCF_000001635.27/genomic.gtf -o 4counts.txt 4.sam
```
save by typing Ctrl x, y, enter

run 

```
sbatch counter.sh
```
</details>
