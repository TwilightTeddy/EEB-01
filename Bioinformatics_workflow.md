## Bioinformatics Pipeline: Trim Galore → STAR → featureCounts  

Below is a breakdown of the core workflow used in our RNA-seq analysis pipeline:

We began with adapter and quality trimming using trimGalore, which removed low-quality bases and discarded reads shorter than 50 bp. This produced cleaned paired-end FASTQ files ready for alignment.

Next, we used STAR to align the trimmed reads to the *Ursus arctos* (brown bear) reference genome (GCA_023065955.2). We retained only uniquely mapped reads, and STAR output sorted BAM files containing the aligned sequences for each sample.

We then used featureCounts to quantify gene expression. It counted how many reads overlapped with exons based on the GTF annotation file, grouping them by 'gene_id' to produce a gene-level count matrix.

For downstream analysis, we used DESeq2 to identify differentially expressed genes. We visualized results through PCA plots and heatmaps, allowing us to observe patterns of upregulated and downregulated genes during hibernation.

Finally, we performed Gene Ontology (GO) and KEGG pathway enrichment analysis to interpret the biological significance of the differentially expressed genes.
