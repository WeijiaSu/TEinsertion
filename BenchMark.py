import pandas as pd
from Bio import SeqIO

pd.set_option("display.max_columns",40)

ref="/data/zhanglab/Weijia_Su/CommonDataSet/TE_full.fa"
ref=list(SeqIO.parse(ref,"fasta"))
names=[rec.id for rec in ref]
TElen=dict(zip([rec.id for rec in ref],[len(str(rec.seq)) for rec in ref]))

def TLDRhelper(filename):
	f=pd.read_table(filename)
	f=f.loc[f["IsSpanRead"]==True]
	f=f.loc[f["Useable"]==True]
	return list(set(f["ReadName"]))


def compareTLDR(TLDRdirPre,TEins,paf):
	TLDRfile=TLDRdirPre+".table.txt"
	f1=pd.read_table(TEins)
	f2=pd.read_table(TLDRfile)
	f2=f2.loc[f2["Family"].isin(names)]
	f2["TElen"]=f2["Family"].apply(lambda x: TElen[x])
	
	f2=f2.loc[f2["Filter"]=="PASS"]
	f2=f2.loc[(f2["StartTE"]<100) & (f2["EndTE"]>f2["TElen"]-100) & (abs(f2["LengthIns"])-abs(f2["TElen"])<=200)]
	print(f1[0:10])
	print(f1.shape)
	print(f2[0:10])
	print(f2.shape)
	print(f2["SpanReads"].sum())	
	l1=set(f1["Readname"])

	TLDRdir=TLDRdirPre+"/"
	
	l2=[]
	for i in list(f2["UUID"]):
		filename=TLDRdir+i+".detail.out"
		l2=l2+TLDRhelper(filename)

	l1=set(l1)
	l2=set(l2)
	print(len(l1))
	print(len(l2))

	U=[u for u in l1 if u in l2]
	l1_u=[i for i in l1 if i not in l2]
	l2_u=[j for j in l2 if j not in l1]

	print(len(U))
	print(len(l1_u))
	print(len(l2_u))
	
	f_paf=pd.read_table(paf,header=None,sep=" ")
	f_paf=f_paf.loc[f_paf[0].isin(l2_u)]
	f_paf=f_paf.groupby([0,5],as_index=False).filter(lambda x:len(x)==1)
	print(f_paf[0:10])
	print(f_paf.shape)

compareTLDR("/data/zhanglab/Weijia_Su/TLDR/Fly/TestPipeline/barcode21.fastq_tldr","/data/zhanglab/Weijia_Su/TestInsertion/test_fl_fullLen_insertion.tsv","/data/zhanglab/Weijia_Su/TestInsertion/barcode21.fastq_TE_full.fa.paf")
