# Tool Name: featureCounts
featureCounts is a highly efficient general-purpose read summarization program that counts mapped reads for genomic features such as genes, exons, promoter, gene bodies, genomic bins and chromosomal locations. It can be used to count both RNA-seq and genomic DNA-seq reads.

## Purpose
In our workflow, featureCounts was used to quantify gene expression by counting how many aligned RNA-seq fragments (from STAR) overlapped with known exons from a reference GTF annotation.

## Function
featureCounts takes aligned sequencing reads (BAM and GTF annotation files) and determines how many of them fall within predefined regions (like genes or exons) in a genome. These exon-level counts were then summarized using the `gene_id` attribute to produce gene-level expression values. We used this gene counts matrix (gene_id and sample columns, etc.) as input for differential expression analysis in DESeq2. 

## Usage

We used featureCounts in a SLURM script to generate a gene-level count matrix from all aligned BAM files. 

```
#!/bin/bash

#slurm script to run FeatureCounts on STAR aligned BAM files
#SBATCH --job-name=runFeatureCounts		# Job name
#SBATCH --partition=128x24				# Partition name
#SBATCH --mail-type=ALL               		# Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=UCSC_ID@ucsc.edu   	# Where to send mail
#SBATCH --time=0-10:00:00 				# Wall clock time limit in Days-Hours:min:seconds
#SBATCH --ntasks=1                 		# Run a single task
#SBATCH --cpus-per-task=8                	# Use 8 threads for STAR
#SBATCH --output=/hb/groups/sip_eeb_01/name/scripts/logs/runfeaturecounts_%j.out    # Standard output and error log
#SBATCH --error=/hb/groups/sip_eeb_01/name/scripts/logs/runfeaturecounts_%j.err     # Standard output and error log
#SBATCH --mem=8G                    # Allocate memory for the job.

featureCounts -p -a data/genome/GCF_023065955.2_UrsArc2.0_genomic.gtf -F 'GTF' -g gene_id -t exon -T 8 -o analysis/3_featurecounts/bear_adipose_rawCounts.txt analysis/2_star/*.bam

```
--p: Tells featureCounts that the input data is paired-end

--a <string>: Provides name of an annotation file. 

--F (isGTFAnnotationFile): Specify the format of the annotation file.

--g < string > (GTF.attrType): Specify the attribute type used to group features (eg. exons) into meta-features (eg. genes) when GTF annotation is provided. Ours is gene_id. 

--t < string > (GTF.featureType): Specify the feature type(s) (thread count). 

--o < string >: Output file which saves the gene count matrix to this file. Located in the analysis/3_featurecounts/ directory.

--analysis/2_star/*.bam: For input files. All BAM files from STAR alignment are stored in the analysis/2_star/ folder.




