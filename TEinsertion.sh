#!/bin/bash

read=$1
readName=$(basename $read);
ref=$2
TE=$3
refName=$(basename $ref);
TEname=$(basename $TE)
minimap2 -x map-ont $ref $read -Y -t 16 | awk '{for(i=1;i<=9;i++) printf $i" "; print ""}' FS='\t' > $readName"_"$refName".paf";
minimap2 -x map-ont $TE $read -Y -t 16 | awk '{for(i=1;i<=9;i++) printf $i" "; print ""}' FS='\t' > $readName"_"$TEname".paf";
python3 /data/zhanglab/Weijia_Su/Git/TEinsertion/TEinsertion_v2.py -Ga $readName"_"$refName".paf" -Ta $readName"_"$TEname".paf" -pName $readName 

