# Tool name
Sequence Read Archive (SRA) is a library of raw DNA and RNA sequencing data.


## Purpose
SRA's purpose is to 


## Function
Fastq data is highly compressed, which is why we can’t just download fastq files directly. To get data from the SRA, we use two tools: prefetch, which downloads .sra files in the raw compressed format used by NCBI and fasterq-dump, which converts .sra files into .fastq files, which are text files of sequencing reads. Each read has a name, DNA sequence, and quality score.

## Usage
Our .sra files are containers for raw sequencing reads
Explain each parameter used

