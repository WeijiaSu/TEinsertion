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

TE="/data/zhanglab/Weijia_Su/CommonDataSet/TE_full.fa"
dm6="/data/zhanglab/Weijia_Su/Genomes/Dro/dm6.fa"
#reads="/data/zhanglab/Weijia_Su/Nanopore_Raw_Data/171107_LW1/Fly_HMS-Beagle-GFP_TEactive_gDNA.fastq"
reads="/data/zhanglab/Weijia_Su/Nanopore_Raw_Data/20220323_fly_gDNA_driver_lines/barcode21.fastq"
outputName="test_fl"

bash /data/zhanglab/Weijia_Su/Git/TEinsertion/TEinsertion_FL.sh $reads $TE $dm6 $outputName
