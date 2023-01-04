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
parser.add_argument("-t","--TE_map")
parser.add_argument("-r","--Genome_map")
parser.add_argument("-flex","--flexibility",default=100)

args=parser.parse_args()

TEmap=args.TE_map
Refmap=args.Genome_map
offset=args.flexibility

def getFullLength(TEmap,Refmap):
	te=pd.read_table(TEmap,header=None,sep=" ")
	te=te.loc[(te[7]<=offset) & (te[8]>=te[6]-offset)]
	te=te[range(9)]
	ge=pd.read_table(Refmap,header=None,sep=" ")
	ge=ge[range(9)]
	ge=ge.loc[ge[0].isin(te[0])]
	te.columns=["ReadName","ReadLen","Readte_s","Readte_e","TEStrand","TEName","TELen","TE_s","TE_e"]
	ge.columns=["ReadName","ReadLen","Readref_s","Readref_e","RefStrand","RefName","RefLen","Ref_s","Ref_e"]
	f=ge.merge(te,on=["ReadName","ReadLen"],how="inner")
	f=f.loc[(f["Readref_s"]<f["Readte_s"]-offset) & (f["Readref_e"]>f["Readte_e"]+offset)]
	print(f.shape)
	print(f[0:10])
	print(len)
getFullLength(TEmap,Refmap)
