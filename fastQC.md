# Tool Name: fastQC
fastQC is a quality control tool for high-throughput sequencing data. It provides summary statistics and visual diagnostics to help assess the quality of raw or processed sequencing reads. FastQC can detect common issues like poor base quality, contamination, overrepresented sequences, bias, and duplicated reads.

## Purpose
In our workflow, fastQC was used to assess the quality of RNA-seq reads before and after trimming (trimGalore). It helped confirm that adapter sequences and low-quality bases were successfully removed by Trim Galore and that the remaining reads met quality standards for alignment and downstream analysis.

## Function  
fastQC performs a series of quality control checks on FASTQ files and reports metrics. 

## Usage  

We used the following command to trim our paired-end reads and run FastQC (see trimgalore.md for more details):

```
trim_galore --paired -q 24 --fastqc \
--fastqc_args "--noextract --nogroup --outdir 1_trim/fastqc" \
--stringency 5 --illumina --length 50 \
-o 1_trim --clip_R1 12 --clip_R2 12 \
[path/to/read1] [path/to/read2]
```

--paired: Indicates input files are paired-end (two reads per fragment). Read 1 (R1) starts from the 5' end of the fragment, while Read 2 (R2) starts from the 3' end.

-q 24: Trims low-quality bases with Phred scores below 24 from both ends of reads.

--fastqc: Runs FastQC automatically on the trimmed reads.

--fastqc_args "--noextract --nogroup --outdir 1_trim/fastqc":

--noextract: Keeps the FastQC results as .zip files (doesn't unzip them).

--nogroup: Disables per-tile aggregation, which can hide poor quality regions in large files.

--outdir 1_trim/fastqc: Stores FastQC reports in the specified folder.

--stringency 5: Requires a minimum of 5 bp overlap between the adapter and the read for trimming to occur.

--illumina: Looks specifically for Illumina adapter sequences.

--length 50: Discards reads shorter than 50 bp after trimming.

-o 1_trim: Outputs all trimmed FASTQ files to the 1_trim directory.

--clip_R1 12: Removes the first 12 bases of Read 1.

--clip_R2 12: Removes the first 12 bases of Read 2.

[path/to/read1] [path/to/read2]: Input raw FASTQ files from NCBI.

```

