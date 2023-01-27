import pandas as pd

pd.set_option("display.max_columns",40)

ref="/data/zhanglab/Weijia_Su/CommonDataSet/TE_full.fa"
TElen=dict(zip([rec.id for rec in ref],[len(str(rec.seq)) for rec in ref]))

def compareTLDR(TLDRfile,TEins):
	f1=pd.read_table(TEins)
	f2=pd.read_table(TLDRfile)
	f2["TElen"]=f2["Family"].apply(lambda x: TElen[x])
	f2=f2.loc[f2["Filter"]=="PASS"]
	#f2=f2.loc[(f2["StartTE"]<100) & (f["EndTE"]>)
	print(f1[0:10])
	print(f1.shape)
	print(f2[0:10])
	print(f2.shape)

compareTLDR("/data/zhanglab/Weijia_Su/TLDR/Fly/TestPipeline/barcode21.fastq_tldr.table.txt","/data/zhanglab/Weijia_Su/TestInsertion/test_fl_fullLen_insertion.tsv")
