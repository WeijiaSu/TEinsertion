import pandas as pd
import os
import os.path
import argparse
import re
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import time
from cigar import Cigar
from collections import Counter
import pysam
import warnings
warnings.filterwarnings('ignore')


pd.set_option("display.max_column",40)

parser=argparse.ArgumentParser()
parser.add_argument("-Ta","--TE_bam")
parser.add_argument("-Ga","--Genome_bam")
parser.add_argument("-pName","--Prefix")
parser.add_argument("-fq","--Rawfastq")
parser.add_argument("-genome","--Genome")
parser.add_argument("-flex","--flexibility",default=100)

args=parser.parse_args()

Ta=args.TE_bam
Ga=args.Genome_bam
pName=args.Prefix
fq=args.Rawfastq
fl=args.flexibility
genome=args.Genome

def FullLegnthTE(TEfile):
	# Read paf of TE and the genome alignment
	TE=pd.read_table(TEfile,header=None,sep=" ")
	# Get the first 9 columns of both files.
	TE=TE[range(0,9)]
	# Re-name the columns
	TE_columns=["Readname","ReadLen","ReadStart_TE","ReadEnd_TE","Strand_TE","TE_Name","TElen","TEstart","TEend"]
	TE.columns=TE_columns
	
	#Select full length
	TE=TE.loc[(TE["TEstart"]<=fl) & (TE["TEend"]>=TE["TElen"]-fl)]
	TE=TE.loc[(TE["ReadStart_TE"]>=100) | (TE["ReadEnd_TE"]<=TE["ReadLen"]-100)]
	left=TE[["Readname","ReadStart_TE"]]
	left["s"]=0
	left=left[["Readname","s","ReadStart_TE"]]
	right=TE[["Readname","ReadEnd_TE","ReadLen"]]
	left.to_csv(pName+".left",header=None,index=None,sep="\t")
	right.to_csv(pName+".right",header=None,index=None,sep="\t")
	print(TE.shape)
	print(TE[0:10])

#FullLegnthTE(Ta)

def getSeq(fastq,Genome):
	seqtk1="seqtk subseq %s %s > %s"%(fastq,pName+".left",pName+".left.fastq")
	seqtk2="seqtk subseq %s %s > %s"%(fastq,pName+".right",pName+".right.fastq")
	os.system(seqtk1)
	os.system(seqtk2)
	minimap="minimap2 -x map-ont %s %s -Y -t 16 | awk '{for(i=1;i<=9;i++) printf $i\" \";print \"\"}' FS='\\t' > %s"%(Genome,pName+".left.fastq",pName+".left.paf")
	os.system(minimap)
	minimap="minimap2 -x map-ont %s %s -Y -t 16 | awk '{for(i=1;i<=9;i++) printf $i\" \";print \"\"}' FS='\\t' > %s"%(Genome,pName+".right.fastq",pName+".right.paf")
	os.system(minimap)
	print(minimap)
#getSeq(fq,genome)

def getJuctions(paf1,paf2):
	left=pd.read_table(paf1,header=None,sep=" ")
	left[10]=left[3]-left[2]
	left=left.drop_duplicates([0,2,3],keep=False).sort_values([0,10],ascending=[True,False])
	left=left.drop_duplicates([0],keep="first")
	left=left.loc[abs(left[1]-left[3])<=500]
	right=pd.read_table(paf2,header=None,sep=" ")
	right[10]=right[3]-right[2]
	right=right.drop_duplicates([0,2,3],keep=False).sort_values([0,10],ascending=[True,False])
	right=right.drop_duplicates([0],keep="first")
	right=right.loc[right[2]<=500]
	left[9]=left[0].apply(lambda x:x.split(":")[0])
	right[9]=right[0].apply(lambda x:x.split(":")[0])
	combined=pd.merge(left,right,on=[9],how="inner")
	combined=combined[(combined["5_x"]==combined["5_y"]) & (combined["4_x"]==combined["4_y"])]
	combined["J1"]=0
	combined["J2"]=0

	combined.loc[combined["4_x"]=="+","J1"]=combined["8_x"]
	combined.loc[combined["4_x"]=="-","J1"]=combined["7_x"]
	combined.loc[combined["4_y"]=="+","J2"]=combined["7_y"]
	combined.loc[combined["4_y"]=="-","J2"]=combined["8_y"]
	combined=combined.loc[abs(combined["J1"]-combined["J2"])<1000]
	print(combined[0:10])
	print(combined.shape)
	print(combined.drop_duplicates([9],keep="first").shape)
	combined.to_csv(pName+".insReads.tsv",header=None,index=None)
getJuctions(pName+".left.paf",pName+".right.paf")
