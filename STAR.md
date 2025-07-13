# Tool Name: STAR
STAR (Spliced Transcripts Alignment to a Reference) is an aligner designed to specifically address many of the challenges of RNA-seq data mapping using a strategy to account for spliced alignments. In other words, STAR splits each read into chunks, finds where those chunks match the genome, trusts the ones that map uniquely, groups nearby matches, and then connects them into a full read alignment — scoring each possibility to pick the best one. STAR is shown to have high accuracy and outperforms other aligners by more than a factor of 50 in mapping speed, but it is memory intensive. 


## Purpose
For every read that STAR aligns, STAR will search for the longest sequence that exactly matches one or more locations on the reference genome. These longest matching sequences are called MMPs. The different parts of the read that are mapped separately are called ‘seeds’. The algorithm achieves this highly efficient mapping by performing a two-step process including seed searching, clustering, stitching, and scoring. 

## Function
STAR indexes genomes, which is a very efficient function as it does not require scanning the entire genome for every read. STAR reads your input genome FASTA file(s), which contain the DNA sequences of chromosomes/scaffolds. Then, it combines these sequences into a single internal representation. In our workflow, we used STAR to align the trimmed reads to the Ursus arctos (brown bear) reference genome (GCA_023065955.2). We retained only uniquely mapped reads, and STAR output sorted BAM files containing the aligned sequences for each sample.

## Usage

We created a SLURM script to handle both downloading and converting STAR data in one batch job. Below is an example of our SLURM script and an explanation of each part:

```
#slurm script to run STAR alignment on trimmed fastq data
#SBATCH --job-name=runSTAR		# Job name
#SBATCH --partition=128x24				# Partition name
#SBATCH --mail-type=ALL               		# Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=UCSC_ID@ucsc.edu   	# Where to send mail
#SBATCH --time=0-10:00:00 				# Wall clock time limit in Days-Hours:min:seconds
#SBATCH --ntasks=1                 		# Run a single task
#SBATCH --cpus-per-task=12                	# Use 12 threads for STAR
#SBATCH --output=/hb/groups/sip_eeb_01/name/scripts/logs/index_genome_star_%j.out    # Standard output and error log
#SBATCH --error=/hb/groups/sip_eeb_01/name/scripts/logs/index_genome_star_%j.err     # Standard output and error log
#SBATCH --mem=32G                    # Allocate memory for the job.
```
Below are the command lines we used to run STAR:
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





