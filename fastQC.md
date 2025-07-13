# Tool Name: fastQC
fastQC is a quality control tool for high-throughput sequencing data. It provides summary statistics and visual diagnostics to help assess the quality of raw or processed sequencing reads. FastQC can detect common issues like poor base quality, contamination, overrepresented sequences, bias, and duplicated reads.

## Purpose
In our workflow, fastQC was used to assess the quality of RNA-seq reads before and after trimming (trimGalore). It helped confirm that adapter sequences and low-quality bases were successfully removed by Trim Galore and that the remaining reads met quality standards for alignment and downstream analysis.

## Function
fastQC performs a series of quality control checks on FASTQ files and reports metrics such as:

## Usage


```

