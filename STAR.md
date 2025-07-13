# Tool Name: STAR
STAR (Spliced Transcripts Alignment to a Reference) is an aligner designed to specifically address many of the challenges of RNA-seq data mapping using a strategy to account for spliced alignments. In other words, STAR splits each read into chunks, finds where those chunks match the genome, trusts the ones that map uniquely, groups nearby matches, and then connects them into a full read alignment — scoring each possibility to pick the best one. STAR is shown to have high accuracy and outperforms other aligners by more than a factor of 50 in mapping speed, but it is memory intensive. 


## Purpose
For every read that STAR aligns, STAR will search for the longest sequence that exactly matches one or more locations on the reference genome. These longest matching sequences are called MMPs. The different parts of the read that are mapped separately are called ‘seeds’. The algorithm achieves this highly efficient mapping by performing a two-step process including seed searching, clustering, stitching, and scoring. 

## Function
STAR indexes genomes, which is a very efficient function as it does not require scanning the entire genome for every read. 

1. STAR reads your input genome FASTA file(s), which contain the DNA sequences of chromosomes/scaffolds. Then, it combines these sequences into a single internal representation and records chromosome names and lengths.

2. STAR builds a suffix array and related data structures that allow for ultra-fast searching of subsequences (e.g., mapping reads to genome locations). It also builds and index of the genome to allow for quick and memory-efficient lookups.

3. Incorporate Transcript Annotations (Optional, but Recommended)
If you provide a --sjdbGTFfile, STAR parses the GTF to find exon-exon splice junctions.
For each annotated junction, STAR extracts a piece of the genome around the junction (--sjdbOverhang, typically Read Length - 1).
These sequences are added to a splice junction database, improving STAR’s ability to align spliced RNA-seq reads that cross exon boundaries.

4. Prepare Output Index Files
STAR saves a set of binary and text files in the --genomeDir, which include:
Genome — compressed version of the genome sequence.
SA and SAindex — suffix array and its index.
chrName.txt — chromosome names.
sjdbList.out.tab — list of splice junctions (if GTF provided).
sjdbInfo.txt — info about splice junctions.
Other internal files. These files are later used in the alignment step when you provide --genomeDir.

## Usage

We created a SLURM script to handle both downloading and converting STAR data in one batch job. Below is an example of our SLURM script and an explanation of each part:

```
STAR \
--runThreadN 12 \
--genomeDir ./2_star/indexed_genome/ \
--sjdbGTFfile ./star/ursus_arctos.gtf \
--sjdbOverhang 100 \
--outFilterMultimapNmax 1 \
--readFilesIn 1_trim/read_1_val.fastq 1_trim/read_2_val.fastq \
--twopassMode Basic \
--outFileNamePrefix analysis/2_star/${sra}_
--outSAMtype BAM SortedByCoordinate
```
--runThreadN: Defines the # of threads to be used for genome generation.

--runMode genomeGenerate option directs STAR to run genome indices generation job. 

-- genomeDir: Specifies path to the directory where the genome indices are stored.

-- genomeFastaFiles: Specifies one or more FASTA files with the genome reference sequences.

-- sjdbGTFfile: Specifies the path to the file with annotated transcripts.

--sjdbOverhang: Specifies the length of the genomic sequence around the annotated junction.

--outFilterMultimapNmax: Read alignments will be output only if the read maps is less or equal than this value, otherwise no alignments will be output

--twopassMode Basic: In the 2-pass mapping job, STAR will map the reads twice. In the 1st pass, the novel junctions will be detected and inserted into the genome indices. In the 2nd pass, all reads will be re-mapped using annotated (from the GTF file) and novel (detected in the 1st pass) junctions.





