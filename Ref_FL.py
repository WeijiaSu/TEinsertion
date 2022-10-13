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

def getFullLength(MapRes,Refmap):
	te=pd.read_table(MapRes,header=None,sep=" ")
	te=te.loc[(te[7]<=offset) & (te[8]>=te[6]-offset)]
	print(te[0:10])
	print(te.shape)
	print(len(set(te[0])))
	ge=pd.read_table(Refmap,header=None,sep=" ")
	ge=ge[range(9)]
	ge=ge.loc[ge[0].isin(te[0])]
	print(ge.shape)
	print(ge[0:10])
	print(len(set(ge[0])))
	
	

	f=ge.merge(te)

getFullLength(TEmap,Refmap)

