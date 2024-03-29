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
pd.set_option("display.max_columns",40)


f1=pd.read_table("/data/zhanglab/Weijia_Su/TestInsertion/test_061223/shmCherry.fastq-TE_full.fa.bam_AligTable.tsv")
f1=f1.loc[f1["Strand"]=="-"]
print(f1[0:20])
print(f1.shape)

f2=pd.read_table("/data/zhanglab/Weijia_Su/TestInsertion/test_061223/shmCherry.fastq-TE_full.fa.bam.paf")
f2=f2.loc[f2["Strand"]=="-"]
print(f2[0:20])
print(f2.shape)

