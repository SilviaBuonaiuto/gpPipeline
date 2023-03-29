#### Required tools and programming language
- bwa (http://bio-bwa.sourceforge.net/bwa.shtml)
- sambamba (https://lomereiter.github.io/sambamba/docs/sambamba-view.html)
- freebayes (https://github.com/freebayes/freebayes)
- samtools/htslib (https://github.com/samtools/htslib)
- bcftools(http://samtools.github.io/bcftools/bcftools.html)
- Vt (https://genome.sph.umich.edu/wiki/Vt)
- vep (https://www.ensembl.org/info/docs/tools/vep/index.html)
- vcftools (https://vcftools.github.io/index.html)
- Python3
- R

#### 0. Index reference genome
index reference sequences in the FASTA format
```
bwa index hg38.p12.fa
```

#### 1. Align reads to reference genome 
```
mkdir bam
bwa mem -t 24 \
-R "@RG\tID:$id\tSM:$id" reference.GRCh38.fa \
sample_R1.fastq.gz sample_R2.fastq.gz | \
samtools view -b - > bam/sample.raw.bam && touch sample.align_ok
```
- sort bed file

```
sambamba sort -t 16 \
-m 32G \
--tmpdir /scratch \
-o sample.raw.sorted.bam \
bam/sample.raw.bam
```
- remove PCR duplicates
```
sambamba markdup -t 8 \
--tmpdir ~/tmp \
--overflow-list-size 500000 \
bam/sample.raw.sorted.bam \
bam/sample.bam
```
#### 2. Variant Calling by sample and chromosome
```
mkdir vcf
freebayes -f reference.GRCh38.fa \
-r $chr \
-g 500 \
-b bam/sample.bam > vcf/sample.$chr.fb.vcf && touch sample.$chr.fb_ok
```
-bgzip and tabix
```
bgzip sample.$chr.fb.vcf
```
```
tabix -p vcf sample.$chr.fb.vcf.gz
```
- filter vcf for quality > 20
```
bcftools filter -i "QUAL>20" \
-O z \
-o vcf/sample.$chr.fb.filtered.vcf.gz \
vcf/sample.$chr.fb.vcf.gz
```
- tabix
- normalize filtered vcf
```
vt normalize -n \
-r GRCh38.fa sample.$chr.filtered.fb.vcf.gz\
> vcf/sample.$chr.fb.norm.vcf && touch sample.$chr_norm_ok
```
- bgzip and tabix
- decompose normalized vcf
```
vt  decompose_blocksub sample.$chr.fb.norm.vcf.gz \
> vcf/sample.$chr.fb.decomp.vcf && touch sample.$chr_decomp_ok
```
-bgzip and tabix

#### 3. Merge all samples
```
bcftools merge -O z \
-0 \
-o testdata/sam/samples.chr22.vcf.gz
-m none \ 
vcf/sample1.$chr.fb.decomp.vcf.gz  vcf/sample2.$chr.fb.decomp.vcf.gz vcf/sample3.$chr.fb.decomp.vcf.gz vcf/sampleX.$chr.fb.decomp.vcf.gz 
```
- bgzip and tabix

#### 4.  Annotate variants using `V.E.P`
merged vcf per chromosome. output is a table

testdata/sam/samples.chr22.vcf.gz -> testdata/sam/samples.chr22.vep.tsv
testdata/con/controls.chr22.vcf -> testdata/con/controls.chr22.vep.tsv  
ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/HGDP/data/

```
vep --af_1kg --af_gnomad --appris --biotype --buffer_size 5000 --check_existing --distance 5000 --fork 4 --polyphen b --pubmed --regulatory --sift b --species homo_sapiens --symbol --tsl --cache --dir_cache /data/biocontainers/vepcache --offline --tab --fields "Uploaded_variation,Location,Allele,Gene,Feature,Feature_type,Consequence,cDNA_position,CDS_position,Protein_position,Amino_acids,Codons,Existing_variation,IMPACT,SYMBOL,STRAND,SIFT,PolyPhen,EXON,AF,AFR_AF,AMR_AF,ASN_AF,EUR_AF,EAS_AF,SAS_AF,AA_AF,EA_AF,gnomAD_AF,gnomAD_AFR_AF,gnomAD_AMR_AF,gnomAD_ASJ_AF,gnomAD_EAS_AF,gnomAD_FIN_AF,gnomAD_NFE_AF,gnomAD_OTH_AF,gnomAD_SAS_AF,MAX_AF,CADD_RAW,CADD_PHRED" --force_overwrite --variant_class -i allsamples.chrx.vcf --plugin CADD,/CADD/whole_genome_SNVs.tsv.gz -o  allsamples.$chr.vep.tsv && touch tabOk/allsamples.$chr.table_ok
```

#### 4.1  Remove Header of file VEP
for i in {1..22} X ; do grep -v "##" ${i}.vep.tsv > ${i}.vep.noheader.tsv; done
for i in {1..22} X; do sed -i '/Uploaded_variatio/s/^#//g' ${i}.vep.noheader.tsv ; done

#### 5. Extract individual's allele count from vcf 
1. Extract counts with vcftools 
```
mkdir counts
cd counts
mkdir $(chr)
vcftools --gzvcf samples.chr22.vcf.gz \
--out $(chr)/$(id).chr22_counts \
--counts \
--indv $(id)
```
2. Change file formatting with [altCounts.py](https://github.com/SilviaBuonaiuto/grepPipeline/tree/master/scr/altCounts.py)
```
-i = path to input file
-id = sample ID
python3 altCounts.py -i /$(chr)/$(id).$(chr)_counts.frq.count -id $(id)
```

#### 6. Create db, add annotations and filter variable sites with [GP.py](https://github.com/SilviaBuonaiuto/gpPipeline/tree/master/scr/GP.py)
*Per chromosome*

```
-db DB                path to samples db
-isa ISA              path to samples vep table
-ic IC                path to controls vep table
-gn GN                path to gene lists file
-p P                  path to pLI score table
-ff FF                define rare frequency treshold
-f F                  define additional rare frequency treshold
-os OS                path to samples output file
-oc OC                path to controls output file
-type TYPE            feature_type (Transcript , intergenic, regularoty)
-r R                  set to false to obtain variants with frequency <= of threshold
-pli PLI              threshold for pLI score
-cadd CADD            treshold for CADD score
-g G                  number of gene lists
-cl CL                list of control samples id
-i I                  number of iterations
-n N                  number of individual to sample
-pathTodirCtrl PATHTODIRCTRL  path to file with control allele counts directory
-chrom CHROM          chromosome name
-ac AC                allele count >= of
-ctgn CTGN            path to control genes file
-gtd GTD              path to genes to discard file
-gt GT                threshold for excluding genes
-sl SL                list of sampes id
-pathTodir PATHTODIR  path to sample allele counts directory
-o O                  path to output file

```
*Command line example*
```
python3 scr/GP.py -db chr22.db -isa testdata/sam/samples.chr22.vep.tsv -ic controls.chr22.vep.tsv -gn testdata/db/all_geneList.tsv -p testdata/db/pLIscore.tsv -ff 0.01 -f 0.05 -os samples.chr22.csq.tsv -oc controls.chr22.csq.tsv -type Transcript -r false -pli 0.7 -cadd 0.5 -g 2 -cl testdata/con/controls_id.txt -i 100 -n 10 -pathTodirCtrl /control/alleleCount/chr22 -chrom chr22 -ac 1 -ctgn controlGenes.chr22.tsv -gtd GenesToDiscardHgdp.chr22.tsv -gt 0.5 -sl testdata/sam/samples_id.txt -pathTodir /samples/counts/chr22 -o samples.chr22.filtered.tsv
```
*Concat all chromosomes -> allSamples.filtered.tsv*

#### 7. Format results with [grep.curation_nolow.R](https://github.com/SilviaBuonaiuto/gpPipeline/tree/master/scr/grep.curation_nolow.R) 
merge all chromosomes samples files
```
Options : 
filtered samples file
plot and tables names prefix

command line example :
Rscript grepPipeline/scr/grep.curation_nolow.R allSamples.filtered.tsv Grep_allsamples
```
