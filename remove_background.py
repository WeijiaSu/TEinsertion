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

Chromosome=["chr2L","chr2R","chr3L","chr3R","chr4","chrX","chrY"]

parser=argparse.ArgumentParser()
parser.add_argument("-I1","--file1",help='first file')
parser.add_argument("-I2","--file2",help="second file")

args=parser.parse_args()

f1_=args.file1
f2_=args.file2

def orgFile(file1):
	f=pd.read_table(file1)
	f["chrome"]=0
	f.loc[f["RName_ref1"]=="0","chrome"]=f["RName_ref2"]
	f.loc[f["RName_ref1"]!="0","chrome"]=f["RName_ref1"]
	f["J"]=0
	f.loc[f["J1"]==0,"J"]=f["J2"].apply(lambda x: round(x,-3))
	f.loc[f["J1"]!=0,"J"]=f["J1"].apply(lambda x: round(x,-3))
	print(f[0:10])
	print(f.shape)
	return f

def removeBG(f1_,f2_):
	f1=orgFile(f1_)
	f2=orgFile(f2_)
	
	print(f1.shape)
	print(f1[0:10])

	print(f2.shape)
	print(f2[0:10])

removeBG(f1_,f2_)
