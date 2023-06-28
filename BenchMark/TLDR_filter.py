import pandas as pd

pd.set_option("display.max_columns",40)

def getFile(name):
	f=pd.read_table(name)
	#f=f.loc[f["Family"]=="HMS-Beagle"]
	f=f.loc[f["Filter"]=="PASS"]
	#f=f.loc[f["StartTE"]<=100 & (f["EndTE"]>=7062-100)]
	print(f.shape)
	print(f[0:10])


for i in ["shwhite","shmCherry","shaubago3"]:
	getFile(i+".fastq_dm6.fa_tldr.table.txt")
