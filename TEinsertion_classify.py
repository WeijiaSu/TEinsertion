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
parser.add_argument("-pName","--Prefix")

args=parser.parse_args()

pName=args.Prefix

def Classfy_TEmap(TE_paf):
	f=pd.read_table(TE_paf,header=None,sep=" ")
	f1=f.loc[(f[7]<=100) & (f[8]>=f[6]-100)]
	f1.to_csv(TE_paf+"_fulllength.tsv",header=None,index=None,sep="\t")
	f2=f.loc[((f[7]<=100) & (f[8]<f[6]-100)) | ((f[7]>100) & (f[8]>=f[6]-100))]
	f2.to_csv(TE_paf+"_OneEnd.tsv",header=None,index=None,sep="\t")
	f3=f.loc[(f[7]>100) & (f[8]<f[6]-100)]
	f3=f3.to_csv(TE_paf+"_NoEnd.tsv",header=None,index=None,sep="\t")
	
	print(f[0:10])
	print(f.shape)


Classfy_TEmap("barcode21.fastq_TE_full.fa.paf")
