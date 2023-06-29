import pandas as pd
from Bio import SeqIO

# set the maximum number of columns displayed in the console to 40
pd.set_option("display.max_columns", 40)

# store the reference file location in the variable "ref"
ref = "/data/zhanglab/Weijia_Su/CommonDataSet/TE_full.fa"
# parse the reference file and store the record IDs in the "names" list
ref = list(SeqIO.parse(ref, "fasta"))
names = [rec.id for rec in ref]
# create a dictionary where the key is the record ID and the value is the length of the sequence
TElen = dict(zip([rec.id for rec in ref], [len(str(rec.seq)) for rec in ref]))

# function to read the table from the TLDR file and filter the usable reads
def TLDRhelper(filename):
	f = pd.read_table(filename)
	# select only rows where "IsSpanRead" is True and "Useable" is True
	f = f.loc[f["IsSpanRead"] == True]
	f = f.loc[f["Useable"] == True]
	return f


# function to compare the TE insertion and TLDR results
def compareTLDR(TLDRfile):
	f2 = pd.read_table(TLDRfile)
	f2 = f2.loc[f2["Family"].isin(names)]
	# add a new column "TElen" to the "f2" dataframe based on the length of the sequence
	f2["TElen"] = f2["Family"].apply(lambda x: TElen[x])
	l=["opus", "copia", "micropia", "mdg3", "Quasimodo", "412"]
	f2=f2.loc[f2["Family"].isin(l)]
	# and "StartTE" less than 100, "EndTE" greater than the sequence length minus 100,
	# and the difference between "LengthIns" and the sequence length being less than or equal to 200
	f2 = f2.loc[f2["Filter"] == "PASS"]
	f2 = f2.loc[(f2["StartTE"] < 100) & (f2["EndTE"] > f2["TElen"] - 100)]
	g=f2.groupby(["Family"],as_index=False).count().sort_values(["SampleReads"],ascending=[False])[["Family","SampleReads"]]
	
	# print the first 10 rows of the "f2" dataframe and its shape
	#print(f2[0:10])
	print(f2.shape)
	#print(f2["SpanReads"].sum())
	print(g[0:10])


for i in ["shwhite","shmCherry","shaubago3"]:
	compareTLDR(i+".fastq_dm6.fa_tldr.table.txt")
