# Tool Name: TrimGalore

## Purpose  
TrimGalore is used for **automated adapter trimming** and **quality filtering** of FASTQ files prior to alignment. It ensures that low-quality or contaminated sequences are removed, improving alignment accuracy and downstream quantification.  
FastQC is integrated to generate quality control reports of the trimmed reads.

## Function  
TrimGalore performs the following operations on paired-end sequencing data:

- Detects and removes Illumina adapter sequences
- Trims low-quality ends from both reads (Phred score < 24)
- Clips the first 12 bases from each read regardless of quality
- Discards reads that are shorter than 50 base pairs after trimming
- Runs FastQC on the trimmed output to evaluate read quality
- Ensures adapter trimming only occurs when a minimum overlap of 5 bp exists between read and adapter

## Usage  
We used the following command to trim our paired-end reads and run FastQC:

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


