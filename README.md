# TE-Genome Analysis Pipeline

This script provides a comprehensive pipeline for analyzing the relationship between transposable elements (TEs) and their respective genomes using mapped read data.

## Dependencies

Ensure that you have the following dependencies installed:

- `pandas`
- `os`
- `argparse`
- `re`
- `Bio` from `biopython`
- `time`
- `cigar`
- `collections`
- `pysam`
- `warnings`
- `bamToPaf` (custom module)

## Features

- Extract mapped reads from BAM files.
- Map those extracted reads to a genome using `minimap2`.
- Convert BAM alignment to PAF format.
- Filter out reads from the TE PAF and Genome PAF.
- Merge the TE and Genome alignments.
- Derive insertions based on the alignments.
- Extract junctions for single and double insertions.
- Combine single and double junction data to produce a final insertion list.

## Usage

```bash
python mean_script.py -t [TE_bam] -n [Name Prefix] -fq [Rawfastq] -g [Genome] -flex [Flexibility]

Arguments
-t, --TE_bam: BAM file with reads mapped to transposon consensus sequences.
-n, --Name: Prefix for output files.
-fq, --Rawfastq: Provide the raw fastq file.
-g, --Genome: Provide the reference genome.
-flex, --flexibility: Set the flexibility parameter (default is 100).

# Workflow Steps
Extract mapped reads from the BAM file using samtools.
Map these reads to the genome.
Convert the resulting BAM alignment to PAF format.
Filter the reads from the TE PAF to ensure that only meaningful reads are used for subsequent analysis.
Similarly, filter the reads from the Genome PAF.
Merge the TE and Genome alignments.
Identify and categorize the insertions into single, double, and multiple.
Extract junctions for single insertions.
Extract junctions for double insertions.
Combine the results of the single and double junction extractions.

# Output Files
The script generates several output files, such as:

.mappedTE.fastq: Extracted mapped reads.
.genome.bam: Reads mapped to genome.
_TE.paf: PAF of TE reads.
_genome.paf: PAF of genome reads.
_merged.tsv: Merged alignments of TE and genome.
.single.tsv, .double.tsv, .multiple.tsv: Categorized insertions.
.single_junction.tsv, .double_junction.tsv: Extracted junctions.
.insertion.tsv: Final insertion list.


