#!/bin/bash

## Modify this job script accordingly and submit with the command:
##    sbatch HPC.sbatch
#SBATCH --nodes=1   # number of nodes
#BATCH --ntasks-per-node=1   # 16 processor core(s) per node
#SBATCH --job-name='Insertion'
#SBATCH --mem=100000
#SBATCH --partition="all"
#SBATCH --mail-user=weijia.su@duke.edu   # email address
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --output="Insertion-%j.out"
#SBATCH --error="Insertion-%j.err"
## LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE


############################################################################################################################
###  This script can be used to detect all full length copies (both referenced and non referenced) ;
###  It needs 4 parameters:
###	 1. The transposon consensus sequences in fasta format (In the example, I used all fly transposons)
###	 2. The genome sequence (for example, dm6.fa)
###  3. The raw reads in fastq format.
###  4. A name for output prefix. (For example: test...)
###  The pipeline will first do preprocessing for the raw reads. Then map them to the TE consensus, select the full length reads, and locate the break poin
###  To run the pipeline, you can change the paprameters (lines 30-33), keep the last line (#35 unchanged)


TE="/data/zhanglab/Weijia_Su/CommonDataSet/TE_full.fa"
dm6="/data/zhanglab/Weijia_Su/Genomes/Dro/dm6.fa"
reads="/data/zhanglab/Weijia_Su/Nanopore_Raw_Data/20220323_fly_gDNA_driver_lines/barcode21.fastq"
outputName="test_fl"

bash /data/zhanglab/Weijia_Su/Git/TEinsertion/TEinsertion_FL.sh $reads $TE $dm6 $outputName
