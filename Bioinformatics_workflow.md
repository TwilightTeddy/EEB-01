#Bioinformatics pipeline:  trimgalore → STAR → featureCounts

adapt and quality trimming with trimgalore (remove <50 bp reads), got FASTQ files from this, mapping transcripts to brown bear genome with STAR <- Map to brown bear genome (GCA_023065955.2)
  Keep only uniquely mapped reads
Convert to BAM
Sort for downstream use,

quantifying gene counts with featureCounts <- Count reads overlapping exons
Assign to genes using GTF

THEN WE GET GENE_COUNTMATRIX

differential gene expression analysis with DESeq2, exploratory analysis with heatmap and PCA, then gene ontology and KEGG enrichment analysis

