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
	#f_te["Q_S"]=f_te["QName"].apply(str)+"_"+f_te["RName"].apply(str)
	bedA=f_te[["QName","QLen"]].drop_duplicates(["QName"],keep="first")
	bedA["s"]=0
	bedA[["QName","s","QLen"]].to_csv(pName+".bedA.bed",header=None,index=None,sep="\t")
	f_te[["QName","QStart","QEnd"]].to_csv(pName+".bedB.bed",header=None,index=None,sep="\t")
	bedtools="bedtools coverage -a %s -b %s > %s"%(pName+".bedA.bed",pName+".bedB.bed",pName+".bed")
	os.system(bedtools)
	
	f_bed=pd.read_table(pName+".bed",header=None)
	f_bed["cutoff"]=f_bed[2]-f_bed[2]*f_bed[6]
	f_bed=f_bed.loc[f_bed["cutoff"]>=fl*2]
	f_te=f_te.loc[f_te["QName"].isin(f_bed[0])]
	f_te.to_csv(pName+"_TE.paf"+".filter.paf",index=None,sep="\t")
	os.remove(pName+".bedA.bed")
	os.remove(pName+".bedB.bed")
	os.remove(pName+".bed")

def filterGenomeReads(Ge_paf):
	f_ge=pd.read_table(Ge_paf)
	f_ge=f_ge.sort_values(["QName","QStart"])
	f_ge_full=f_ge.loc[(f_ge["QStart"]<fl) & (f_ge["QEnd"]>f_ge["QLen"]-fl)]
	f_ge=f_ge.loc[~f_ge["QName"].isin(f_ge_full["QName"])]
	f_ge=f_ge.loc[f_ge["RName"].isin(Chromosome)]
	f_ge.to_csv(pName+"_genome.paf"+".filter.paf",index=None,sep="\t")


def combineAlignment(TE_paf,Ge_paf):
	f_te=pd.read_table(TE_paf)
	f_ge=pd.read_table(Ge_paf)
	f=f_ge.merge(f_te,on=["QName"],how="inner")
	f=f.loc[(f["QEnd_x"]<f["QStart_y"]+fl) | (f["QStart_x"]>f["QEnd_y"]-fl)]
	#f=f.loc[(abs(f["QEnd_x"]-f["QStart_y"])<=5*fl) | (abs(f["QStart_x"]-f["QEnd_y"])<=5*fl)]
	f["overlap"]=False
	overlap1=((f["QStart_x"]<=f["QStart_y"]) & (f["QEnd_x"]>=f["QEnd_y"]))
	overlap2=((f["QStart_y"]<=f["QStart_x"]) & (f["QEnd_y"]>=f["QStart_x"]))
	overlap3=((f["QEnd_y"]>=f["QEnd_x"]) & (f["QStart_y"]<=f["QEnd_x"]-fl))
	f.loc[overlap1,"overlap"]="overlap1"
	f.loc[overlap2,"overlap"]="overlap2"
	f=f.loc[f["overlap"]==False]
	f=f.loc[(abs(f["QEnd_x"]-f["QStart_y"])<=5*fl) | (abs(f["QStart_x"]-f["QEnd_y"])<=5*fl)]
	f=f.drop(["overlap"],axis=1)
	f.to_csv(pName+"_merged.tsv",index=None,sep="\t")

def getInsertion(Filename):
	f=pd.read_table(Filename)
	f=f.sort_values(["QName","QStart_x","QEnd_x"])
	f=f.drop_duplicates(["QName","QStart_x","QEnd_x"],keep="first")
	f=f.drop_duplicates(["QName","QStart_x"],keep="last")
	f=f.drop_duplicates(["QName","QEnd_x"],keep="first")
	f1=f.groupby(["QName","RName_y","QStart_y","QEnd_y","RStart_y","REnd_y"],as_index=False).filter(lambda x: len(x)==1)
	f2=f.groupby(["QName","RName_y","QStart_y","QEnd_y","RStart_y","REnd_y"],as_index=False).filter(lambda x: len(x)==2)
	f3=f.loc[~(f["QName"].isin(f1["QName"])) & ~(f["QName"].isin(f2["QName"]))]
	f1.to_csv(pName+".single.tsv",index=None,sep="\t")
	f2.to_csv(pName+".double.tsv",index=None,sep="\t")
	f3.to_csv(pName+".multiple.tsv",index=None,sep="\t")


def getSingle(single):
	f=pd.read_table(single)
	f["ds1"]=abs(f["QEnd_x"]-f["QStart_y"])
	f["ds2"]=abs(f["QStart_x"]-f["QEnd_y"])
	f["J1"]=0
	f["J2"]=0
	f.loc[(f["ds1"]<f["ds2"])&(f["Strand_x"]=="+"),"J1"]=f["REnd_x"]
	f.loc[(f["ds1"]<f["ds2"])&(f["Strand_x"]=="-"),"J1"]=f["RStart_x"]
	f.loc[(f["ds1"]>=f["ds2"])&(f["Strand_x"]=="+"),"J2"]=f["RStart_x"]
	f.loc[(f["ds1"]>=f["ds2"])&(f["Strand_x"]=="-"),"J2"]=f["REnd_x"]
	f["flanking"]="single"
	#f=f.drop(["ds1","ds2"],axis=1)
	#f["QStart_ref1"]=0
	#f["QEnd_ref1"]=0
	#f["QStart_TE"]=f["RStart_y"]
	#f["QEnd_TE"]=f["QEnd_y"]
	#f["QStart_ref2"]=0
	#f["QEnd_ref2"]=0
	#f["RName_ref1"]=0
	#f["RLen_ref"]=0
	#f["RStart_ref1"]=0
	#f["REnd_ref1"]=0
	#f["Strand_ref1"]=0
	#f["RName_TE"]=f["RName_y"]
	#f["RLen_TE"]=f["RLen_y"]
	#f["RStart_TE"]=f["RStart_y"]
	#f["REnd_TE"]=f["REnd_y"]
	#f["Strand_TE"]=f["Strand_y"]
	#f["RName_ref2"]=0
	#f["RStart_ref2"]=0
	#f["REnd_ref2"]=0
	#f["Strand_ref2"]=0

	#f.loc[(f["ds1"]<f["ds2"])&(f["Strand_x"]=="+"),"QStart_ref1"]=f["QStart_x"]
	#f.loc[(f["ds1"]<f["ds2"])&(f["Strand_x"]=="+"),"QEnd_ref1"]=f["QEnd_x"]
	#f.loc[(f["ds1"]<f["ds2"])&(f["Strand_x"]=="+"),"QStart_ref2"]="NA"
	#f.loc[(f["ds1"]<f["ds2"])&(f["Strand_x"]=="+"),"QEnd_ref2"]="NA"
	#f.loc[(f["ds1"]<f["ds2"])&(f["Strand_x"]=="+"),"RName_ref1"]=f["RName_x"]
	#f.loc[(f["ds1"]<f["ds2"])&(f["Strand_x"]=="+"),"RLen_ref"]=f["RLen_x"]
	#f.loc[(f["ds1"]<f["ds2"])&(f["Strand_x"]=="+"),"RStart_ref1"]=f["RStart_x"]
	#f.loc[(f["ds1"]<f["ds2"])&(f["Strand_x"]=="+"),"REnd_ref1"]=f["REnd_x"]
	#f.loc[(f["ds1"]<f["ds2"])&(f["Strand_x"]=="+"),"REnd_ref1"]=f["REnd_x"]
	#


	#
	f.to_csv(pName+".single_junction.tsv",index=None,sep="\t")
	#
	

def getDouble(double):
	f=pd.read_table(double)
	f=f.drop(["QLen_y"],axis=1)
	f=f.sort_values(["QName","QStart_x","QEnd_x"])
	f_s=f.drop_duplicates(["QName"],keep="first")
	f_e=f.drop_duplicates(["QName"],keep="last")
	f_new=f_s.merge(f_e,on=["QName","QStart_y","QEnd_y","QLen_x"])
	f_new=f_new.loc[f_new["RName_x_x"]==f_new["RName_x_y"]]
	f_new=f_new.loc[f_new["Strand_x_x"]==f_new["Strand_x_y"]]
	
	f_new["ds1"]=abs(f_new["QEnd_x_x"]-f_new["QStart_y"])
	f_new["ds2"]=abs(f_new["QStart_x_y"]-f_new["QEnd_y"])
	f_new=f_new.loc[(f_new["ds1"]<=fl*5) & (f_new["ds2"]<=fl*5)]

	f_new=f_new[["QName","QLen_x","QStart_x_x","QEnd_x_x","QStart_y","QEnd_y","QStart_x_y","QEnd_x_y","RName_x_x","RLen_x_x","RStart_x_x","REnd_x_x","Strand_x_x","RName_y_x","RLen_y_x","RStart_y_x","REnd_y_x","Strand_y_x","RName_x_y","RStart_x_y","REnd_x_y","Strand_x_y"]]
	columns=["QName","QLen","QStart_ref1","QEnd_ref1","QStart_TE","QEnd_TE","QStart_ref2","QEnd_ref2","RName_ref1","RLen_ref","RStart_ref1","REnd_ref1","Strand_ref1","RName_TE","RLen_TE","RStart_TE","REnd_TE","Strand_TE","RName_ref2","RStart_ref2","REnd_ref2","Strand_ref2"]
	f_new.columns=columns

	f_new["J1"]=0
	f_new["J2"]=0

	f_new.loc[f_new["Strand_ref1"]=="+","J1"]=f_new["REnd_ref1"]
	f_new.loc[f_new["Strand_ref1"]=="-","J1"]=f_new["RStart_ref1"]
	f_new.loc[f_new["Strand_ref2"]=="+","J2"]=f_new["RStart_ref2"]
	f_new.loc[f_new["Strand_ref2"]=="-","J2"]=f_new["REnd_ref2"]
	
	f_new=f_new.loc[abs(f_new["J1"]-f_new["J2"])<=fl*5]
	f_new.to_csv(pName+".double_junction.tsv",index=None,sep="\t")
	f_new["flanking"]="double"


def	classifySingle(SingleFile):
	f=pd.read_table(SingleFile)
	f=f.loc[(f["RStart_y"]<=fl) | (f["REnd_y"]>=f["RLen_y"]-fl)]
	f=f.groupby(["QName"],as_index=False).filter(lambda x : len(x)==1)	
	g=f.groupby(["RName_y"],as_index=False).count().sort_values(["QName"],ascending=False)[["RName_y","QName"]]
	
	#l=["opus", "copia", "micropia", "mdg3", "Quasimodo", "412"]
	#g=g.loc[g["RName_y"].isin(l)]
	print(g[0:10])
	print(f.shape)
	print(f.drop_duplicates(["QName"],keep="first").shape)


def classifyDouble(doubleFile):
	print(pName+".double_junction.tsv")
	f=pd.read_table(doubleFile)
	f=f.loc[(f["RStart_TE"]<=fl) & (f["REnd_TE"]>=f["RLen_TE"]-fl)]
	f=f.groupby(["QName"],as_index=False).filter(lambda x : len(x)==1)
	g=f.groupby(["RName_TE"],as_index=False).count().sort_values(["QName"],ascending=False)[["RName_TE","QName"]]
	l=["opus", "copia", "micropia", "mdg3", "Quasimodo", "412"]
	f=f.loc[f["RName_TE"].isin(l)]
	f=f.sort_values(["RName_TE"])
	#g=g.loc[g["RName_TE"].isin(l)]
	#print(g[0:10])
	print(f.shape)
	print(f)
	f.to_csv(pName+".reads.tsv",index=None,sep="\t")

#getMappedReads(Ta)
#MapToGenome()
#convertToPaf(Ta,pName+"_TE")
#convertToPaf(pName+"_genome.bam",pName+"_genome")
#filterTEreads(pName+"_TE.paf")
#filterGenomeReads(pName+"_genome.paf")
#combineAlignment(pName+"_TE.paf"+".filter.paf",pName+"_genome.paf"+".filter.paf")
#getInsertion(pName+"_merged.tsv")
#getSingle(pName+".single.tsv")
#getDouble(pName+".double.tsv")
#classifySingle(pName+".single_junction.tsv")
classifyDouble(pName+".double_junction.tsv")
