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
parser.add_argument("-flex","--flexibility",default=100)

args=parser.parse_args()

Ta=args.TE_bam
Ga=args.Genome_bam
pName=args.Prefix
fq=args.Rawfastq
fl=args.flexibility


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

def getSeq(fastq):
	#seqtk1="seqtk subseq %s %s > %s"%(fastq,pName+".left",pName+".left.fastq")
	#seqtk2="seqtk subseq %s %s > %s"%(fastq,pName+".right",pName+".right.fastq")
	#os.system(seqtk1)
	#os.system(seqtk2)
	minimap="minimap2 -x map-ont $ref $reads -Y -t 16 | awk '{for(i=1;i<=9;i++) printf $i" ";}' FS='\\t'"
	print(minimap)
getSeq(fq)