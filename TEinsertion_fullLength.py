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
parser.add_argument("-pName","--Prefix")
parser.add_argument("-fq","--Rawfastq")
parser.add_argument("-genome","--Genome")
parser.add_argument("-buf","--bufferLen",default=100)

args=parser.parse_args()

Ta=args.TE_bam
pName=args.Prefix
fq=args.Rawfastq
fl=args.bufferLen
genome=args.Genome

def FullLegnthTE(TE):
	
	#Select full length
	TE_Name=list(TE["TE_Name"])[0]
	# Filtering TE to only include rows where TEstart is less or equal to fl
	# and TEend is greater or equal to TElen minus fl
	TE=TE.loc[(TE["TEstart"]<=fl) & (TE["TEend"]>=TE["TElen"]-fl)]
	
	# Further filtering of TE to only include rows where ReadStart_TE is greater or equal to 100
	# or ReadEnd_TE is less or equal to ReadLen minus 100
	TE=TE.loc[(TE["ReadStart_TE"]>=100) | (TE["ReadEnd_TE"]<=TE["ReadLen"]-100)]
	
	
	if TE.shape[0]>0:
		# Prepare left and right dataframes for CSV export
		left=TE[["Readname","ReadStart_TE"]]
		left["s"]=0
		left=left[["Readname","s","ReadStart_TE"]]
		# Export left and right dataframes to CSV
		right=TE[["Readname","ReadEnd_TE","ReadLen"]]
		left.to_csv(pName+".left",header=None,index=None,sep="\t")
		right.to_csv(pName+".right",header=None,index=None,sep="\t")
		print("Full length TE: %s"%(TE.shape[0]))
	else:
		print("No full length %s"%(TE_Name))
	return TE

def getSeq(fastq,Genome):
	# Define the seqtk command line for the left.fastq and right.fastq output
	seqtk1="seqtk subseq %s %s > %s"%(fastq,pName+".left",pName+".left.fastq")
	seqtk2="seqtk subseq %s %s > %s"%(fastq,pName+".right",pName+".right.fastq")
	os.system(seqtk1)
	os.system(seqtk2)

	#Map the letf and right fastq files to the genome using minimap2 and save the results to paf files
	minimap="minimap2 -x map-ont %s %s -Y -t 16 | awk '{for(i=1;i<=9;i++) printf $i\" \";print \"\"}' FS='\\t' > %s"%(Genome,pName+".left.fastq",pName+".left.paf")
	os.system(minimap)
	minimap="minimap2 -x map-ont %s %s -Y -t 16 | awk '{for(i=1;i<=9;i++) printf $i\" \";print \"\"}' FS='\\t' > %s"%(Genome,pName+".right.fastq",pName+".right.paf")
	os.system(minimap)

def getJuctions(paf1,paf2):
	left=pd.read_table(paf1,header=None,sep=" ")
	left[10]=left[3]-left[2]
	left=left.drop_duplicates([0,2,3],keep=False).sort_values([0,10],ascending=[True,False])
	left=left.drop_duplicates([0],keep="first")
	left=left.loc[abs(left[1]-left[3])<=500]
	right=pd.read_table(paf2,header=None,sep=" ")
	right[10]=right[3]-right[2]
	right=right.drop_duplicates([0,2,3],keep=False).sort_values([0,10],ascending=[True,False])
	right=right.drop_duplicates([0],keep="first")
	right=right.loc[right[2]<=500]
	left[9]=left[0].apply(lambda x:x.split(":")[0])
	right[9]=right[0].apply(lambda x:x.split(":")[0])
	combined=pd.merge(left,right,on=[9],how="inner")
	combined=combined[(combined["5_x"]==combined["5_y"]) & (combined["4_x"]==combined["4_y"])]
	combined["J1"]=0
	combined["J2"]=0

	combined.loc[combined["4_x"]=="+","J1"]=combined["8_x"]
	combined.loc[combined["4_x"]=="-","J1"]=combined["7_x"]
	combined.loc[combined["4_y"]=="+","J2"]=combined["7_y"]
	combined.loc[combined["4_y"]=="-","J2"]=combined["8_y"]
	combined=combined.loc[abs(combined["J1"]-combined["J2"])<1000]
	if combined.shape[0]>0:
		print("flanking %s"%(combined.shape[0]))
		combined.to_csv(pName+".insReads.tsv",header=None,index=None,sep="\t")
	return combined

def Cluster(combined_result):
	f=pd.read_table(combined_result,header=None)
	bedInput=f[[16,21,22,9]]
	bedInput.loc[bedInput[21]<=bedInput[22],"J1"]=bedInput[21]
	bedInput.loc[bedInput[21]<=bedInput[22],"J2"]=bedInput[22]
	bedInput.loc[bedInput[21]>=bedInput[22],"J1"]=bedInput[22]
	bedInput.loc[bedInput[21]>=bedInput[22],"J2"]=bedInput[21]
	bedInput=bedInput[[16,"J1","J2",9]]
	bedInput["J1"]=bedInput["J1"].apply(int)
	bedInput["J2"]=bedInput["J2"].apply(int)
	bedInput=bedInput.sort_values([16,"J1","J2",9])
	bedInput.to_csv(pName+"bedInput.tsv",header=None,index=None,sep="\t")
	bedtools="bedtools cluster -i %s -d 1000 > %s"%(pName+"bedInput.tsv",pName+"cluster.tsv")
	os.system(bedtools)


def Assign(Cluster,combined,fl_TE):
	f=pd.read_table(combined,header=None)
	c=pd.read_table(Cluster,header=None)
	f.loc[f[21]<=f[22],"J1"]=f[21]
	f.loc[f[21]<=f[22],"J2"]=f[22]
	f.loc[f[21]>=f[22],"J1"]=f[22]
	f.loc[f[21]>=f[22],"J2"]=f[21]
	f["J1"]=f["J1"].apply(int)
	f["J2"]=f["J2"].apply(int)
	c["UID"]=c[3]+"—"+c[0]+"_"+c[1].apply(str)+"_"+c[2].apply(str)
	f["UID"]=f[9]+"—"+f[16]+"_"+f["J1"].apply(str)+"_"+f["J2"].apply(str)
	d=dict(zip(c["UID"],c[4]))
	f["cluster"]=f["UID"].apply(lambda x: d[x])
	f=f.sort_values(["cluster"]).drop_duplicates([9],keep="first")
	f=f[[9,16,"J1","J2","cluster"]]
	f.columns=["Readname","REFname","REFstart","REFend","cluster"]
	fl_TE=fl_TE.drop_duplicates(["Readname"],keep="first")
	fl_TE=fl_TE[["Readname","TE_Name"]]
	res=f.merge(fl_TE,on=["Readname"],how="inner")
	res=res.sort_values(["cluster"])
	return res



def Main(TEmappingDf):
	full_TE=FullLegnthTE(TEmappingDf)
	TEres=pd.DataFrame(columns=["Readname","REFname","REFstart","REFend","cluster","TE_Name"])
	if full_TE.shape[0]>0:
		getSeq(fq,genome)
		com=getJuctions(pName+".left.paf",pName+".right.paf")
		if com.shape[0]>0:
			Cluster(pName+".insReads.tsv")
			TEres=Assign(pName+"cluster.tsv",pName+".insReads.tsv",full_TE)
			#os.remove(pName+".left")
			#os.remove(pName+".right")
			#os.remove(pName+".left.fastq")
			#os.remove(pName+".right.fastq")
			#os.remove(pName+".left.paf")
			#os.remove(pName+".right.paf")
			#os.remove(pName+".insReads.tsv")
			#os.remove(pName+"bedInput.tsv")
			#os.remove(pName+"cluster.tsv")
	return TEres


final_res=pd.DataFrame(columns=["Readname","REFname","REFstart","REFend","cluster","TE_Name"])
TEmap=pd.read_table(Ta,header=None,sep=" ")
TEmap=TEmap[range(0,9)]
TEmap_columns=["Readname","ReadLen","ReadStart_TE","ReadEnd_TE","Strand_TE","TE_Name","TElen","TEstart","TEend"]
TEmap.columns=TEmap_columns
TElist=list(TEmap.drop_duplicates(["TE_Name"],keep="first")["TE_Name"])
for i in TElist:
	print("################ %s %s ####################" %(TElist.index(i)+1,i))
	sub=TEmap.loc[TEmap["TE_Name"]==i]
	print(sub.shape)
	res=Main(sub)
	print(res[0:10])
	print("################ %s %s ####################" %(i, res.shape[0]))
	final_res=final_res.append(res)
	print("current output")
	print(final_res.shape)
final_res.to_csv(pName+"_fullLen_insertion.tsv",index=None,sep="\t")

