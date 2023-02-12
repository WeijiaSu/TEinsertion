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
def compareTLDR(TLDRdirPre, TEins, paf):
    # store the TLDR file location in the "TLDRfile" variable
    TLDRfile = TLDRdirPre + ".table.txt"
    # read the TE insertion file into the "f1" dataframe
    f1 = pd.read_table(TEins)
    # read the TLDR file into the "f2" dataframe and filter it to only include rows with "Family" in "names"
    f2 = pd.read_table(TLDRfile)
    f2 = f2.loc[f2["Family"].isin(names)]
    # add a new column "TElen" to the "f2" dataframe based on the length of the sequence
    f2["TElen"] = f2["Family"].apply(lambda x: TElen[x])

    # further filter the "f2" dataframe to only include rows with "Filter" value of "PASS",
    # and "StartTE" less than 100, "EndTE" greater than the sequence length minus 100,
    # and the difference between "LengthIns" and the sequence length being less than or equal to 200
    f2 = f2.loc[f2["Filter"] == "PASS"]
    f2 = f2.loc[(f2["StartTE"] < 100) & (f2["EndTE"] > f2["TElen"] - 100) & (abs(f2["LengthIns"]) - abs(f2["TElen"]) <= 200)]

    # print the first 10 rows of the "f1" dataframe and its shape
    print(f1[0:10])
    print(f1.shape)
    # print the first 10 rows of the "f2" dataframe and its shape
    print(f2[0:10])
    print(f2.shape)



	print(f2["SpanReads"].sum())	
	l1=set(f1["Readname"])

	TLDRdir=TLDRdirPre+"/"
	
	l2=[]
	for i in list(f2["UUID"]):
		filename=TLDRdir+i+".detail.out"
		l2.append(TLDRhelper(filename))
	f_t=l2[0]
	for fi in l2[1:]:
		f_t=f_t.append(fi)
	f_t.to_csv(TLDRdirPre+".tldrReads.txt",index=None,sep="\t")

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
