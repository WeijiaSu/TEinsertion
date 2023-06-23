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

Classfy_TEmap("barcode21.fastq_TE_full.fa.paf")
