## Bioinformatics Pipeline:  trimgalore → STAR → featureCounts

Here is an explanation of a key pathway in our workflow: our bioinformatics pipeline. 

We started off by performing adaptation and quality trimming with trimgalore (which removed <50 bp reads), and got FASTQ files as an output. 

Then, we proceeded by mapping transcripts to the brown bear (Ursus arctos) genome with STAR (GCA_023065955.2), keeping only uniquely mapped reads which were converted to BAM files representing the aligned sequences (sorted for downstream use).

Furthemore, we quantified gene counts with featureCounts (reads overlapping exons and assigns to genes using GTF), which outputted a gene count matrix. 

We then conducted differential gene expression analysis with DESeq2 by generating heatmaps, which we used to determine the genes which were upregulated or downregulated during hibernation as well as PCA. 

Finally, we did gene ontology and KEGG enrichment analysis. 
