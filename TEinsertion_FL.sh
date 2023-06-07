#!/bin/bash

reads=$1
TE=$2
ref=$3
outPre=$4

readName=$(basename $reads);
porechop -i $reads --extra_end_trim 0 --discard_middle -o $readName".pre.fastq"
refName=$(basename $ref);
TEname=$(basename $TE)
minimap2 -x map-ont $TE $readName".pre.fastq" -Y -t 16 | awk '{for(i=1;i<=9;i++) printf $i" "; print ""}' FS='\t' > $readName"_"$TEname".paf";
#python3 /data/zhanglab/Weijia_Su/Git/TEinsertion/TEinsertion_fullLength.py -Ta $readName"_"$TEname".paf" -fq $readName".pre.fastq" -genome $ref -pName $outPre

