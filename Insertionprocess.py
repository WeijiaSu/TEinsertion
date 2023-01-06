

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
 parser.add_argument("-flex","--flexibility",default=100)
 
 args=parser.parse_args()
 
 Ta=args.TE_bam
 Ga=args.Genome_bam
 pName=args.Prefix
 fl=args.flexibility




def genomeLocation(outfile):
	f=pd.read_table(outfile)
	f["coor1"]=-1
	f["coor2"]=-1
	f.loc[f["Junc_2"]==-1,"coor1"]=f["Junc_1"].apply(lambda x: round(x,-2))
	f.loc[f["Junc_1"]==-1,"coor1"]=f["Junc_2"].apply(lambda x: round(x,-2))
	f.loc[(f["Junc_1"]!=-1)&(f["Junc_2"]!=-1) & (f["Junc_1"]<=f["Junc_2"]), "coor1"]=f["Junc_1"].apply(lambda x: round(x,-2))
	f.loc[(f["Junc_1"]!=-1)&(f["Junc_2"]!=-1) & (f["Junc_1"]>f["Junc_2"]), "coor1"]=f["Junc_2"].apply(lambda x: round(x,-2))
	f.loc[(f["Junc_1"]!=-1)&(f["Junc_2"]!=-1) & (f["Junc_1"]<=f["Junc_2"]), "coor2"]=f["Junc_2"].apply(lambda x: round(x,-2))
	f.loc[(f["Junc_1"]!=-1)& (f["Junc_2"]!=-1)& (f["Junc_1"]>=f["Junc_2"]), "coor2"]=f["Junc_1"].apply(lambda x: round(x,-2))
	f.loc[f["left_refName"]!="-1","refName"]=f["left_refName"]
	f.loc[f["left_refName"]=="-1","refName"]=f["right_refName"]
	f=f.loc[f["TEconf"]=="FL_TE"]
	f=f.groupby(["refName","coor1"],as_index=False).count()
	print(f.shape)
	print(f[0:10])
genomeLocation(pName+"_InsReads.tsv"
