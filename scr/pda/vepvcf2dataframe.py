import pandas as pd
import argparse

def parse_vep_vcf():
    """
    Parse a vep-annotated vcf file, returns a tsv file with one row per encoded transcript
    """
    parser = argparse.ArgumentParser()
    ###~~~ input files
    parser.add_argument("-i", help="path to input vep vcf file ",type=str,required=True)
    ###~~~ output file
    parser.add_argument("-o", help="path to output tsv file",type=str, required= True)
    args = parser.parse_args()

    VEP_columns = ['Allele','Consequence','IMPACT','SYMBOL','Gene','Feature_type','Feature','BIOTYPE','EXON','INTRON','HGVSc','HGVSp','cDNA_position','CDS_position','Protein_position','Amino_acids','Codons','Existing_variation','DISTANCE','STRAND','FLAGS','VARIANT_CLASS','SYMBOL_SOURCE','HGNC_ID','TSL','APPRIS','SIFT','PolyPhen','AFR_AF','AMR_AF','EAS_AF','EUR_AF','SAS_AF','gnomAD_AF','gnomAD_AFR_AF','gnomAD_AMR_AF','gnomAD_ASJ_AF','gnomAD_EAS_AF','gnomAD_FIN_AF','gnomAD_NFE_AF','gnomAD_OTH_AF','gnomAD_SAS_AF','CLIN_SIG','SOMATIC','PHENO','PUBMED','MOTIF_NAME','MOTIF_POS','HIGH_INF_POS','MOTIF_SCORE_CHANGE','CADD_PHRED','CADD_RAW']


    #read vcf from file
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
    df_final.rename(columns={'ALT': 'Allele'}).to_csv(args.o, sep = "\t", index = False)


if __name__ == "__main__":
	parse_vep_vcf()