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
