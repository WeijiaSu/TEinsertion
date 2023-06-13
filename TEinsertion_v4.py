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

def readAlignment(TE_paf,Ge_paf):
	f_te=pd.read_table(TE_paf)
	f_te=f_te.sort_values(["QName","QStart"])
	f_te_full=f_te.loc[(f_te["QStart"]<fl) & (f_te["QEnd"]>f_te["QLen"]-fl)]
	f_te=f_te.loc[~f_te["QName"].isin(f_te_full["QName"])]
	print(f_te.shape)
	print(f_te[0:10])
	f_ge=pd.read_table(Ge_paf)
	f_ge=f_ge.sort_values(["QName","QStart"])
	f_ge_full=f_ge.loc[(f_ge["QStart"]<fl) & (f_ge["QEnd"]>f_ge["QLen"]-fl)]
	f_ge=f_ge.loc[~f_ge["QName"].isin(f_ge_full["QName"])]
	print(f_ge.shape)
	print(f_ge[0:10])

#getMappedReads(Ta)
#MapToGenome()
#convertToPaf(Ta,pName+"_TE")
#convertToPaf(pName+"_genome.bam",pName+"_genome")

readAlignment(pName+"_TE.paf",pName+"_genome.paf")

