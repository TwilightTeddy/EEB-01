# Tool Name: SRA
Sequence Read Archive (SRA) is a library of raw DNA and RNA sequencing data.


## Purpose
SRA's purpose is to store and provide public access to raw sequencing data generated by high-throughput sequencing technologies.


## Function
Fastq data is highly compressed, which is why we cannot download fastq files directly. To get data from the SRA, we use two tools: prefetch, which downloads .sra files in the raw compressed format used by NCBI and fasterq-dump, which converts .sra files into .fastq files, which are text files of sequencing reads. Each read has a name, DNA sequence, and quality score. 

## Usage

We created a SLURM script to handle both downloading and converting SRA data in one batch job. Below is an example of our SLURM script and an explanation of each part:

```
#!/bin/bash

#SBATCH --job-name=getSRA    			# Job name
#SBATCH --partition=128x24				# Partition name
#SBATCH --mail-type=ALL               		# Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=UCSC_ID@ucsc.edu   	# Where to send mail
#SBATCH --time=0-05:00:00 				# Wall clock time limit in Days-Hours:min:seconds
#SBATCH --ntasks=1                    		# Run a single task
#SBATCH --cpus-per-task=4                  	# Use 4 threads for fasterq-dump
#SBATCH --output=scripts/logs/fasterq-dump_%j.out    # Standard output and error log
#SBATCH --error=scripts/logs/fasterq-dump_%j.err     # Standard output and error log
#SBATCH --mem=8G                    		# Allocate memory for the job.
#SBATCH --array=1-11					# array job

```



