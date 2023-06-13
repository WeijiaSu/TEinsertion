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
from bamToPaf import bamConverter 
warnings.filterwarnings('ignore')


pd.set_option("display.max_column",40)

parser=argparse.ArgumentParser()
parser.add_argument("-t","--TE_bam",help='bam file with reads mapped to transposon consensus sequences')
parser.add_argument("-n","--Name",help="Prefix")
parser.add_argument("-fq","--Rawfastq")
parser.add_argument("-g","--Genome")
parser.add_argument("-flex","--flexibility",default=100)

args=parser.parse_args()

Ta=args.TE_bam
pName=args.Name
fq=args.Rawfastq
fl=args.flexibility
genome=args.Genome


def getMappedReads(bam):
	samtools="samtools view -b -F 4 %s | samtools fastq > %s"%(bam,pName+".mappedTE.fastq")
	os.system(samtools)

def MapToGenome():
	minimap2="minimap2 -ax map-ont %s %s -Y -t 16 | samtools view -b | samtools sort > %s"%(genome,pName+".mappedTE.fastq",pName+"_genome.bam")
	os.system(minimap2)

def convertToPaf(bamfile,name):
	bamConverter().ConverAlignment(bamfile,name)

def filterTEreads(TE_paf):
	f_te=pd.read_table(TE_paf)
	f_te["Q_S"]=f_te["QName"].apply(str)+"_"+f_te["RName"].apply(str)
	bedA=f_te[["Q_S","QLen"]].drop_duplicates(["Q_S"],keep="first")
	bedA["s"]=0
	bedA[["Q_S","s","QLen"]].to_csv(pName+".bedA.bed",header=None,index=None,sep="\t")
	f_te[["Q_S","QStart","QEnd"]].to_csv(pName+".bedB.bed",header=None,index=None,sep="\t")
	bedtools="bedtools coverage -a %s -b %s > %s"%(pName+".bedA.bed",pName+".bedB.bed",pName+".bed")
	os.system(bedtools)
	
	f_bed=pd.read_table(pName+".bed",header=None)
	f_bed["cutoff"]=f_bed[2]-f_bed[2]*f_bed[6]
	f_bed=f_bed.loc[f_bed["cutoff"]>=fl*2]
	f_te=f_te.loc[f_te["Q_S"].isin(f_bed[0])]
	f_te=f_te.drop(["Q_S"],axis=1)
	f_te.to_csv(pName+"_TE.paf"+".filter.paf",index=None,sep="\t")
	os.remove(pName+".bedA.bed")
	os.remove(pName+".bedB.bed")
	os.remove(pName+".bed")

def filterGenomeReads(Ge_paf):
	f_ge=pd.read_table(Ge_paf)
	f_ge=f_ge.sort_values(["QName","QStart"])
	f_ge_full=f_ge.loc[(f_ge["QStart"]<fl) & (f_ge["QEnd"]>f_ge["QLen"]-fl)]
	f_ge=f_ge.loc[~f_ge["QName"].isin(f_ge_full["QName"])]
	f_ge.to_csv(pName+"_genome.paf"+".filter.paf",index=None,sep="\t")


def combineAlignment(TE_paf,Ge_paf):
	f_te=pd.read_table(TE_paf)
	f_ge=pd.read_table(Ge_paf)
	print(f_te.shape)
	print(f_te[0:10])
	print(f_ge[0:10])
	f=f_ge.merge(f_te,on=["QName"],how="inner")
	print(f.shape)
	print(f[0:10])
	f=f.loc[(abs(f["QEnd_x"]-f["QStart_y"])<=5*fl) | (abs(f["QStart_x"]-f["QEnd_y"])<=5*fl)]
	print(f.shape)
	print(f[0:10])

#getMappedReads(Ta)
#MapToGenome()
#convertToPaf(Ta,pName+"_TE")
#convertToPaf(pName+"_genome.bam",pName+"_genome")
#filterTEreads(pName+"_TE.paf")
filterGenomeReads(pName+"_genome.paf")
combineAlignment(pName+"_TE.paf"+".filter.paf",pName+"_genome.paf"+".filter.paf")
