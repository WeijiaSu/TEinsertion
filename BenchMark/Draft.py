import pandas as pd

pd.set_option("display.max_columns",40)

def getFile(name):
	f=pd.read_table(name)
	f=f.loc[f["Family"]=="HMS-Beagle"]
	f=f.loc[f["Filter"]=="PASS"]
	f=f.loc[f["StartTE"]<=100 & (f["EndTE"]>=7062-100)]
	print(f.shape)
	#print(f[0:10])


for i in ["white","aub-white","aub-PolQ","aub-XRCC1"]:
	getFile(i+"_combine.fastq.20X.fastq-dm6.fa.bam_tldr.table.txt")
