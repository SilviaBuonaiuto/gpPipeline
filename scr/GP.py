import sqlite3
import pandas as pd
import re, sys, argparse, random, glob, subprocess

def main():
	parser = argparse.ArgumentParser()
	###~~~~ new sql db name
	parser.add_argument("-db", help="path to samples db",type=str)
	###~~~ input files
	parser.add_argument("-isa", help="path to samples vep table",type=str)
	parser.add_argument("-ic", help="path to controls vep table",type=str)
	###~~~~ annotations tables
	parser.add_argument("-gn", help="path to gene lists file")
	parser.add_argument("-p", help="path to pLI score table",type=str)
	###~~~ frequencies options
	parser.add_argument("-ff", help="define rare frequency treshold",type=float)
	parser.add_argument("-f", help="define additional rare frequency treshold",type=float)
	###~~~ output file
	parser.add_argument("-os", help="path to samples output file",type=str)
	parser.add_argument("-oc", help="path to controls output file",type=str)
	###~~~ filter option
	parser.add_argument("-type", help=" feature_type (Transcript , intergenic, regularoty) ", type=str)
	parser.add_argument("-r", help="set to false to obtain variants with frequency <= of threshold", type=str)
	parser.add_argument("-pli", help="threshold for  pLI score ", type=float)
	parser.add_argument("-cadd", help=" treshold for CADD score ",type=float)
	parser.add_argument("-g", help=" number of gene lists ",type=float)
	###~~~ genotypes and second filter options
	parser.add_argument("-cl", help="list of control samples id ",type=str)
	parser.add_argument("-i", help="number of iterations ",type=int)
	parser.add_argument("-n", help="number of individual to sample ",type=int)
	parser.add_argument("-pathTodirCtrl", help="path to file with control allele counts directory",type=str)
	parser.add_argument("-chrom", help="chromosome name  ",type=str)
	parser.add_argument("-ac", help=" allele count >= of   ",type=int)
	parser.add_argument("-ctgn", help="path to control genes file",type=str)
	parser.add_argument("-gtd", help="path to genes to discard file",type=str)
	parser.add_argument("-gt", help="threshold for excluding genes ",type=float)
	parser.add_argument("-sl", help="list of sampes id ",type=str)
	parser.add_argument("-pathTodir", help="path to sample allele counts directory  ",type=str)
	###~~~ final output
	parser.add_argument("-o", help="path to output file",type=str)
	args = parser.parse_args()
	### create new sql db
	conn = sqlite3.connect(args.db)
	c = conn.cursor()
	c.execute("DROP TABLE IF EXISTS myTableSamples;")
	c.execute("DROP TABLE IF EXISTS myTableControls;")
	c.execute("DROP TABLE IF EXISTS pLItable;")
	c.execute("DROP TABLE IF EXISTS geneList;")
	###### open pLI table and create new table inside db
	df_pli=pd.read_table(args.p, sep="\t")
	df_pli.columns = df_pli.columns.str.strip()
	df_pli.to_sql("pLItable", conn)
	###### open geneList table and create new table inside db
	df_genes=pd.read_table(args.gn, sep="\t")
	df_genes.columns = df_genes.columns.str.strip()
	df_genes['sumGenes']= df_genes.loc[:, df_genes.columns != df_genes.columns[0]].sum(axis=1)
	df_genes.to_sql("geneList", conn)

	###### open samples vep table with pandas and create table inside db
	dfs = pd.read_table(args.isa , sep="\t", index_col= "Uploaded_variation").fillna(0)
	dfsAf=dfs[dfs.columns[dfs.columns.to_series().str.contains('AF')]].replace("-", 0.0)
	dfs.update(dfsAf)
	dfsCadd=dfs[dfs.columns[dfs.columns.to_series().str.contains('CADD')]].replace("-", 0.0)
	dfs.update(dfsCadd)
	dfs.columns = dfs.columns.str.strip()
	dfs.to_sql("myTableS", conn)
	
	###### open controls vep table with pandas and create table inside db
	dfc = pd.read_table(args.ic , sep="\t", index_col= "Uploaded_variation").fillna(0)
	dfcAf=dfc[dfc.columns[dfc.columns.to_series().str.contains('AF')]].replace("-", 0.0)
	dfc.update(dfcAf)
	dfcCadd=dfc[dfc.columns[dfc.columns.to_series().str.contains('CADD')]].replace("-", 0.0)
	dfc.update(dfcCadd)
	dfc.columns = dfc.columns.str.strip()
	dfc.to_sql("myTableC", conn)

	#### create new tables with pLI scores (for samples and controls)
	dict_tables = {'myTableS': 'firstjoinS', 'myTableC': 'firstjoinC'}
	for key, value in dict_tables.items():
		createTable = f'CREATE TABLE {value} AS SELECT * FROM {key} WHERE 1=2 ;'
		alterTable = f'ALTER TABLE {value} ADD pLIscore real ;'
		insert = f'INSERT INTO {value} SELECT DISTINCT {key}.*, pLItable.pLI FROM {key} LEFT JOIN pLItable ON {key}.Feature = pLItable.transcript ;'
		update = f'UPDATE {value} SET pLIscore=0 WHERE pLIscore is null ;'
		c.execute(createTable)
		c.execute(alterTable)
		c.execute(insert)
		c.execute(update)
	conn.commit()

	#### create tables with gene lists info (for samples and controls)
	dict_tables2 = {'firstjoinS' : 'genesjoinS', 'firstjoinC' : 'genesjoinC'}
	for key2, value2 in dict_tables2.items():
		create = f'CREATE TABLE {value2} AS SELECT * FROM {key2} WHERE 1=2;'
		c.execute(create)
		gcolumns = len(df_genes.columns)
		for i in range(1,gcolumns):
			x = f'ALTER TABLE {value2} ADD {df_genes.columns[i]} real;'
			c.execute(x)

	conn.commit()

	my_list = []
	for i in range(1, gcolumns):
		y = f'geneList.{df_genes.columns[i]}'
		my_list.append(y)

	required_lists = ', '.join(my_list)
	for key2, value2 in dict_tables2.items():
		string = f'INSERT INTO {value2} SELECT DISTINCT {key2}.*, {required_lists} FROM {key2} LEFT JOIN geneList ON {key2}.Gene = geneList.ensID;'
		c.execute(string)
		conn.commit()

	gcolumns = len(df_genes.columns)
	for key2, value2 in dict_tables2.items():
		for i in range(1,gcolumns):
			z = f'UPDATE {value2} SET {df_genes.columns[i]}=0 WHERE {df_genes.columns[i]} is null;'
			c.execute(z)
	
	conn.commit()

	
	##### create new table with variants id (in the form chromosome:start:end:alt_allele)
	dict_tables3 = {'genesjoinS' : 'indexTableS', 'genesjoinC' : 'indexTableC'}
	for key3, value3 in dict_tables3.items():
		createT = f'CREATE TABLE {value3} AS SELECT * FROM {key3} WHERE 1=2;'
		alterT = f'ALTER TABLE {value3} ADD index_x text;'
		insertT = f'INSERT INTO {value3} SELECT DISTINCT *, Location ||\':/\'|| Allele FROM {key3};'
		c.execute(createT)
		c.execute(alterT)
		c.execute(insertT)

	conn.commit()

	### Add columns with values for rare frequency threshold
	thr = args.ff
	c.execute("DROP TABLE IF EXISTS freqTable1S;")
	c.execute("DROP TABLE IF EXISTS freqTable1C;")
	dict_tables4 = {'indexTableS' : 'freqTable1S', 'indexTableC' : 'freqTable1C'}
	for key4, value4 in dict_tables4.items():
		createF1 = f'CREATE TABLE {value4} AS SELECT * FROM {key4} WHERE 1=2;'
		alterF1 = f'ALTER TABLE {value4} ADD rare1 real ;'
		c.execute(createF1)
		c.execute(alterF1)
	conn.commit()

	fcolumns = len(dfsAf.columns)
	freq_list = []
	for i in range(0, fcolumns-1):
		y = f'{dfsAf.columns[i]}'
		freq_list.append(y)

	nob_freq = ' =0 AND '.join(freq_list)
	values_freq = ' <=? AND '.join(freq_list)
	for key4, value4 in dict_tables4.items():
		freq_string = f'INSERT INTO {value4} SELECT *, CASE WHEN {nob_freq} =0 THEN "NOB" WHEN {values_freq} <=? THEN "true" ELSE "false" END FROM {key4};'
		c.execute(freq_string, (thr,thr,thr,thr,thr,thr,thr,thr,thr,thr,thr,))
	conn.commit()

	### Add columns with values for additional rare frequency threshold
	thr2 = args.f
	c.execute("DROP TABLE IF EXISTS freqTable2S;")
	c.execute("DROP TABLE IF EXISTS freqTable2C;")
	dict_tables5 = {'freqTable1S' : 'freqTable2S', 'freqTable1C' : 'freqTable2C'}
	for key5, value5 in dict_tables5.items():
		createF2 = f'CREATE TABLE {value5} AS SELECT * FROM {key5} WHERE 1=2;'
		alterF2 = f'ALTER TABLE {value5} ADD rare2 real ;'
		c.execute(createF2)
		c.execute(alterF2)
	conn.commit()

	nob_freq = ' =0 AND '.join(freq_list)
	values_freq = ' <=? AND '.join(freq_list)
	for key5, value5 in dict_tables5.items():
		freq_string = f'INSERT INTO {value5} SELECT *, CASE WHEN {nob_freq} =0 THEN "NOB" WHEN {values_freq} <=? THEN "true" ELSE "false" END FROM {key5};'
		c.execute(freq_string, (thr2,thr2,thr2,thr2,thr2,thr2,thr2,thr2,thr2,thr2,thr2,))
	conn.commit()

	for key5, value5 in dict_tables5.items():
		for i in range(0,fcolumns-1):
			x = f'ALTER TABLE {value5} DROP COLUMN {dfsAf.columns[i]} ;'
			c.execute(x)
	conn.commit()

	### Create new table with CADD score percentiles
	c.execute("DROP TABLE IF EXISTS finalTableS;")
	c.execute("DROP TABLE IF EXISTS finalTableC;")
	dict_tables6 = {'freqTable2S' : 'finalTableS', 'freqTable2C' : 'finalTableC'}
	for key6, value6 in dict_tables6.items():
		createFi = f'CREATE TABLE {value6} AS SELECT * FROM {key6} WHERE 1=2;'
		alterFi = f'ALTER TABLE {value6} ADD caddPercent real ;'
		insertCadd = f'INSERT INTO {value6} SELECT * , ROUND(PERCENT_RANK() OVER (ORDER BY CADD_RAW),3) FROM {key6};'
		c.execute(createFi)
		c.execute(alterFi)
		c.execute(insertCadd)
	conn.commit()

	queryS = "SELECT * FROM finalTableS;"
	queryC = "SELECT * FROM finalTableC;"
	df_filtSS = pd.read_sql_query(queryS, conn).fillna(0)
	df_filtSC = pd.read_sql_query(queryC, conn).fillna(0)
	df_filtSS.to_csv(args.os, sep="\t", index=False)
	df_filtSC.to_csv(args.oc, sep="\t", index=False)

	typ = args.type ; rareThresh = args.r  ; pliscore = args.pli ; caddscore = args.cadd ; numgene = args.g
	final_dict = {'finalTableS' : 'filteredS' , 'finalTableC' : 'filteredC'}
	for final_key, final_value in final_dict.items():
		createFilt = f'CREATE TABLE {final_value} AS SELECT * FROM {final_key} WHERE IMPACT != "MODIFIER" AND Feature_type = ? AND rare2 != ? AND (pLIscore >= ? AND caddPercent >= ? OR sumGenes >= ?);'
		c.execute(createFilt, (typ,rareThresh, pliscore, caddscore, numgene))
	conn.commit()

	extractFiltS = "SELECT * FROM filteredS;"
	extractFiltC = "SELECT * FROM filteredC;"
	df_filtS = pd.read_sql_query(extractFiltS, conn).fillna(0)
	df_filtCtr = pd.read_sql_query(extractFiltC, conn).fillna(0)

	##~~ repeat args.i times on args.n samples from controls
	listControl = [line.rstrip('\n') for line in open(args.cl)]
	cycle=0
	while cycle < args.i:
		cycle+=1
		##~~ choose a random sample form controls of size args.n
		randomSample=random.sample(listControl, args.n)
		##~~ integrate with samples ID and csqAlele count
		control_filtered_allsamples=pd.DataFrame()
		for ss in randomSample:
			tmpdf=df_filtCtr
			tmpSamp=pd.read_table('%s/%s.%s_counts.tsv' %( args.pathTodirCtrl, ss, args.chrom) )
			if "chr" in tmpSamp.loc[0].key:
				tmpSamp.key = tmpSamp.key.str.lstrip("chr")
			tmpSamp = tmpSamp[ (tmpSamp['ALTcount']>=args.ac )]
			tmpSamp["sample"] = ss
			if "chr" in tmpdf.loc[0].index_x:
				tmpdf.index_x = tmpdf.index_x.str.lstrip("chr")
			df=tmpdf.reset_index(drop=True).merge(tmpSamp, left_on='index_x', right_on='key').drop("key",axis=1)
			control_filtered_allsamples=pd.concat([control_filtered_allsamples, df])
		tmpg=control_filtered_allsamples[[ 'SYMBOL', 'ID']].drop_duplicates().groupby('SYMBOL').count().transform(lambda x: x /float(args.n) )
		if cycle==1:  genesPerSample=tmpg  # create genesPerSample at first cycle
		else: genesPerSample=genesPerSample.join(tmpg, on='SYMBOL', how='outer', lsuffix='_genesPerSample', rsuffix='_tmpg').fillna(0)
		
	##~~ make GrandMean over args.i and discard genes that on average shows up in args.gt individuals over args.i iterations	
	genesPerSample['GrandMean']=genesPerSample.sum(axis=1 )/float(args.i)
	genesToDiscardControl=genesPerSample[genesPerSample['GrandMean']> float(args.gt)]
	genesPerSample.to_csv(args.ctgn, sep='\t', index=False)
	genesToDiscardControl.to_csv(args.gtd, sep='\t', index=False)

	cols = []
	count = 1
	for column in genesPerSample.columns:
		if column == 'ID_genesPerSample':
			cols.append(f'Cycle_{count}')
			count+=1
			continue
		cols.append(column)

	genesPerSample.columns = cols

	cols2 = []
	count = 1
	for column in genesPerSample.columns:
		if column == 'ID_tmpg':
			cols2.append(f'Cycle_n{count}')
			count+=1
			continue
		cols2.append(column)

	genesPerSample.columns = cols2
	#### load table genePerSample in grep db
	#genesPerSample = pd.read_table(args.ctgn, sep="\t")
	#genesPerSample.columns = genesPerSample.columns.str.strip()
	genesPerSample.to_sql("genesMean", conn, if_exists="replace")

	listSamples = [line.rstrip('\n') for line in open(args.sl)]
	ssall=pd.DataFrame()
	for ss in listSamples:
		tmpgrep=df_filtS
		tmpSamp=pd.read_table('%s/%s.%s_counts.tsv' %(args.pathTodir, ss, args.chrom) )
		if "chr" in tmpSamp.loc[0].key:
			tmpSamp.key = tmpSamp.key.str.lstrip("chr")
		if "chr" in tmpgrep.loc[0].index_x:
			tmpgrep.index_x = tmpgrep.index_x.str.lstrip("chr")
		tmpSamp = tmpSamp[ (tmpSamp['ALTcount']>= args.ac )]
		df=tmpgrep.reset_index(drop=True).merge(tmpSamp, left_on='index_x', right_on='key').drop("key",axis=1)
		ssall=pd.concat([ssall, df])
		#df_filtS.to_csv("filtropdaslq" , sep="\t", index=False)
		

	ssall.to_sql("grepFilter", conn,if_exists="replace")
	c.execute("DROP TABLE IF EXISTS noCtrlGenes;")
	c.execute("CREATE TABLE noCtrlGenes AS SELECT grepFilter.*, CASE WHEN genesMean.GrandMean IS NULL THEN 0 ELSE genesMean.GrandMean END as GrandMean FROM grepFilter LEFT JOIN genesMean ON grepFilter.SYMBOL = genesMean.SYMBOL;")
	conn.commit()
	query = "SELECT * FROM noCtrlGenes;"
	df_end = pd.read_sql_query(query,conn)
	df_end.to_csv(args.o, sep = "\t", na_rep= "NA", index = False)

if __name__ == "__main__":
	main()