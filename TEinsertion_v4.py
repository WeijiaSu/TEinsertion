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
parser.add_argument("-t","--TE_bam",description='bam file with reads mapped to transposon consensus sequences')
Ta=args.TE_bam
TaName=str(Ta.split("/")[-1])
parser.add_argument("-n","--Name",description="Prefix",default=TaName)
parser.add_argument("-fq","--Rawfastq")
parser.add_argument("-genome","--Genome")
parser.add_argument("-flex","--flexibility",default=100)

args=parser.parse_args()

Ta=args.TE_bam
pName=args.Prefix
fq=args.Rawfastq
fl=args.flexibility
genome=args.Genome

bamConverter().ConverAlignment(Ta)
