import sqlite3
import pandas as pd
import re, sys, argparse
import numpy as np

def main():
	parser = argparse.ArgumentParser()
	###~~~~ new sql db name
	parser.add_argument("-db", help="path to input vep table file ",type=str,required=True)
	###~~~ input files
	parser.add_argument("-i", help="path to input vep table file ",type=str,required=True)
	###~~~~ databases
	parser.add_argument("-g", help="path to gene file list ",required=True)
	parser.add_argument("-pli", help="path to table of pLI score  ",type=str, required= True)
	#parser.add_argument('-fathmmcoding', help="path to fathmm file",type=str, required= True)
	#parser.add_argument('-fathmmnc', help="path to fathmm coding file",type=str, required= True)
	###~~~ output file
	parser.add_argument("-o", help="path to output file  ",type=str, required= True)
	args = parser.parse_args()


	conn = sqlite3.connect(args.db)  ### create new sql db
	c = conn.cursor()
	c.execute("DROP TABLE IF EXISTS myTable;")
	c.execute("DROP TABLE IF EXISTS pLItable;")
	c.execute("DROP TABLE IF EXISTS geneList;")
	c.execute("DROP TABLE IF EXISTS fatmTab;")
	c.execute("DROP TABLE IF EXISTS fatmNCtab;")
	###### open vep table with pandas and create table inside grep.db
	VEP_columns = ['Allele','Consequence','IMPACT','SYMBOL','Gene','Feature_type','Feature','BIOTYPE','EXON','INTRON','HGVSc','HGVSp','cDNA_position','CDS_position','Protein_position','Amino_acids','Codons','Existing_variation','DISTANCE','STRAND','FLAGS','VARIANT_CLASS','SYMBOL_SOURCE','HGNC_ID','TSL','APPRIS','SIFT','PolyPhen','AFR_AF','AMR_AF','EAS_AF','EUR_AF','SAS_AF','gnomAD_AF','gnomAD_AFR_AF','gnomAD_AMR_AF','gnomAD_ASJ_AF','gnomAD_EAS_AF','gnomAD_FIN_AF','gnomAD_NFE_AF','gnomAD_OTH_AF','gnomAD_SAS_AF','CLIN_SIG','SOMATIC','PHENO','PUBMED','MOTIF_NAME','MOTIF_POS','HIGH_INF_POS','MOTIF_SCORE_CHANGE','CADD_PHRED','CADD_RAW']
	df = pd.read_csv(args.i,comment='#',sep='\t',header=None)
	df = df.rename(columns ={0:'CHROM',1:'POS',2:'ID',3:'REF',4:'ALT',5:'QUAL',6:'FILTER',7:'INFO'})
	df['Uploaded_variation'] = df['CHROM']+'_'+df['POS'].astype(str)+'_'+df['REF']+'/'+df['ALT']
	#split VEP CSQ annotation field according to encoded transcript
	df['INFOa']= (df['INFO'].str.split('CSQ=')).apply(lambda x: x[0])
	df['INFOb']= (df['INFO'].str.split('CSQ=')).apply(lambda x: x[1])
	df['INFOb_list'] = df['INFOb'].str.split(',')
	#get one row per transcript
	df_explode = df.explode('INFOb_list').reset_index(drop=True)
	#split again to get single VEP annotation fields
	df_info = df_explode['INFOb_list'].str.split('|').apply(pd.Series)
	df_info.set_axis(VEP_columns, axis=1, inplace=True)
	df_concat = pd.concat([df_explode,df_info],axis=1)
	df2 = df_concat[['Uploaded_variation','CHROM','POS','ID','REF','ALT']+VEP_columns].replace('', '-')
	df2['Location'] = df2['CHROM']+':'+df2['POS'].astype(str)
	df_final = df2[['Uploaded_variation' ,'Location' ,'ALT' ,'Gene' ,'Feature' ,'Feature_type' , 'Consequence' , 'cDNA_position', 'CDS_position', 'Protein_position','Amino_acids' ,'Codons' ,'Existing_variation' , 'IMPACT' , 'SYMBOL', 'STRAND' , 'SIFT' , 'PolyPhen' , 'EXON', 'AFR_AF' , 'AMR_AF' , 'EAS_AF','EUR_AF', 'SAS_AF', 'gnomAD_AF', 'gnomAD_AFR_AF', 'gnomAD_AMR_AF', 'gnomAD_ASJ_AF', 'gnomAD_EAS_AF', 'gnomAD_FIN_AF', 'gnomAD_NFE_AF', 'gnomAD_OTH_AF', 'gnomAD_SAS_AF', 'CADD_RAW', 'CADD_PHRED']]
	df_final.rename(columns={'ALT': 'Allele'})
	#df = pd.read_table(args.i, sep="\t", index_col= "Uploaded_variation").fillna(0)
	df_final[['AFR_AF', 'AMR_AF', 'EAS_AF', 'EUR_AF', 'SAS_AF', 'gnomAD_AF', 'gnomAD_AFR_AF', 'gnomAD_AMR_AF', 'gnomAD_ASJ_AF', 'gnomAD_EAS_AF', 'gnomAD_FIN_AF', 'gnomAD_NFE_AF', 'gnomAD_OTH_AF', 'gnomAD_SAS_AF','CADD_RAW', 'CADD_PHRED']] = df[['AFR_AF', 'AMR_AF', 'EAS_AF', 'EUR_AF', 'SAS_AF', 'gnomAD_AF', 'gnomAD_AFR_AF', 'gnomAD_AMR_AF', 'gnomAD_ASJ_AF', 'gnomAD_EAS_AF', 'gnomAD_FIN_AF', 'gnomAD_NFE_AF', 'gnomAD_OTH_AF', 'gnomAD_SAS_AF', 'CADD_RAW', 'CADD_PHRED']].replace('-', 0.0)
	df_final.columns = df.columns.str.strip()
	df_final.to_sql("myTable", conn)
	###### open pLI table and create new table inside grep.db
	df=pd.read_table(args.pli, sep="\t")
	df.columns = df.columns.str.strip()
	df.to_sql("pLItable", conn)
	###### open geneList table and create new table inside grep.db
	df=pd.read_table(args.g, sep="\t")
	df.columns = df.columns.str.strip()
	df.to_sql("geneList", conn)
		
	###### join myTable and pLItab on common column creating new table

	c.execute('CREATE TABLE firstjoin (Uploaded_variation text,Location text,Allele text,Gene text,Feature text,Feature_type text,Consequence text,cDNA_position integer,CDS_position integer,Protein_position integer,Amino_acids text,Codons text,Existing_variation text,IMPACT text,SYMBOL text,STRAND text,SIFT real,PolyPhen real,EXON integer,AFR_AF real,AMR_AF real,EAS_AF real,EUR_AF real,SAS_AF real,gnomAD_AF real,gnomAD_AFR_AF real,gnomAD_AMR_AF real,gnomAD_ASJ_AF real,gnomAD_EAS_AF real,gnomAD_FIN_AF real,gnomAD_NFE_AF real,gnomAD_OTH_AF real,gnomAD_SAS_AF real,CADD_RAW real, CADD_PHRED real, pLIscore real);')

	c.execute('INSERT INTO firstjoin SELECT myTable.*, pLItable.pLI FROM myTable LEFT JOIN pLItable ON myTable.Feature = pLItable.transcript;')
	conn.commit()
	c.execute('UPDATE firstjoin SET pLIscore=0 WHERE pLIscore is null');
	conn.commit()

	##### join genelist

	c.execute('CREATE TABLE genesjoin (Uploaded_variation text,Location text,Allele text,Gene text,Feature text,Feature_type text,Consequence text,cDNA_position integer,CDS_position integer,Protein_position integer,Amino_acids text,Codons text,Existing_variation text,IMPACT text,SYMBOL text,STRAND text,SIFT real,PolyPhen real,EXON integer,AFR_AF real,AMR_AF real,EAS_AF real,EUR_AF real,SAS_AF real,gnomAD_AF real,gnomAD_AFR_AF real,gnomAD_AMR_AF real,gnomAD_ASJ_AF real,gnomAD_EAS_AF real,gnomAD_FIN_AF real,gnomAD_NFE_AF real,gnomAD_OTH_AF real,gnomAD_SAS_AF real,CADD_RAW real, CADD_PHRED real, pLIscore real, EmbryoDev real, DDD real, Lethal real, Essential real, Misc real, DDR real, Candidate real);')

	c.execute('INSERT INTO genesjoin SELECT firstjoin.*, geneList.EmbryoDev, geneList.DDD, geneList.Lethal, geneList.Essential, geneList.Misc, geneList.DDR, geneList.Candidate FROM firstjoin LEFT JOIN geneList ON firstjoin.Gene = geneList.ensID;')
	conn.commit()

	###### add index_x column

	c.execute('CREATE TABLE indexTable (Uploaded_variation text,Location text,Allele text,Gene text,Feature text,Feature_type text,Consequence text,cDNA_position integer,CDS_position integer,Protein_position integer,Amino_acids text,Codons text,Existing_variation text,IMPACT text,SYMBOL text,STRAND text,SIFT real,PolyPhen real,EXON integer,AFR_AF real,AMR_AF real,EAS_AF real,EUR_AF real,SAS_AF real,gnomAD_AF real,gnomAD_AFR_AF real,gnomAD_AMR_AF real,gnomAD_ASJ_AF real,gnomAD_EAS_AF real,gnomAD_FIN_AF real,gnomAD_NFE_AF real,gnomAD_OTH_AF real,gnomAD_SAS_AF real,CADD_RAW real, CADD_PHRED real, pLIscore real, EmbryoDev real, DDD real, Lethal real, Essential real, Misc real, DDR real, Candidate real, index_x text);')

	c.execute("INSERT INTO indexTable SELECT *, Location ||\':/\'|| Allele FROM genesjoin;")

	conn.commit()

	query = "SELECT * FROM indexTable;" 
	df_final = pd.read_sql_query(query, conn).fillna(0)
	df_final.to_csv(args.o, sep="\t", index=False)

if __name__ == "__main__":
	main()

