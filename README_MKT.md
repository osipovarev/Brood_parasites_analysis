## 1. McDonald-Kreitman analysis

Described as performed for cowbirds (Molothrus ater).
The analysis was similar for all analysed species.
Steps 1.1 - 1.5 were performed on the Harvard cluster

All supplemetary scripts can be found in [here](https://github.com/osipovarev/brood_parasites) or [here](https://github.com/osipovarev/general_file_tools)

### 1.1. Prepare genomic data
```
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/012/460/135/GCA_012460135.2_BPBGC_Mater_1.1/GCA_012460135.2_BPBGC_Mater_1.1_genomic.fna.gz

zcat $fna | cut -d" " -f1 | sed "s/\.[1,2]//g" > $db.fa
```


### 1.2. Call variants with snpArcher

#### prepare sample file
```
db=molAte
GENOME=$BP/$db/Genome/$db.fa

write_samples_OK.py -f $FASTQS -r $GENOME -n Molothrus_ater -a GCA_012460135.2 > workflow/config/samples.molAte.csv
```
#### prepare workflow/config/config.molAte.yaml

#### run snpArcher
```
conda activate snparcher
snakemake --cores 1 --profile ../profiles/slurm --configfile config/config.molAte.yaml
```
#### filter final VCF
```
VCF=${db}_full_raw.vcf.gz

bcftools view -f .,PASS $VCF -a -U -O u | bcftools view -v snps -m2 -M2 -e 'F_MISSING > 0.75 | ref="N" | ALT="." | TYPE~"indel"' -O z -o passed.$VCF
```



### 1.3. Prepare NCBI annotation
```
NCBIFTP=https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/012/460/135/GCF_012460135.2_BPBGC_Mater_1.1/
GTF=GCF_012460135.2_BPBGC_Mater_1.1_genomic.gtf
REPORT=GCF_012460135.2_BPBGC_Mater_1.1_assembly_report.txt
FEATURE=GCF_012460135.2_BPBGC_Mater_1.1_feature_table.txt

wget $NCBIFTP/$GTF.gz
wget $NCBIFTP/$REPORT
wget $NCBIFTP/$FEATURE.gz

gunzip $GTF.gz
gunzip $FEATURE.gz

awk -F"\t" '{print $15"\t"$11}' $FEATURE  | awk 'NF==2{print}' |s -u  | grep -v product > ncbi_isoforms.tsv

renameToHLscaffolds.py -d <(grep -v ^# $REPORT | awk '{print $7","$5}' | sed 's/\.2$//' | sed 's/\.1$//') -c 1 -a $GTF | sed 's/^na/JAAQZI010000148/' > $db.ncbi_anno.gtf
```



### 1.4. Add outgroups to the ingroup VCF

[Call variants from genome alignment recipe](https://github.com/lh3/minimap2/blob/master/cookbook.md#genome-aln))
[minimap2 source code and installation](https://github.com/lh3/minimap2/tree/master)(


#### Align agePho/ageTri (red-winged and tricolored blackbirds) to molAte brown-headed cowbird) with minimap2
```
BP=/n/holylfs05/LABS/informatics/Lab/project-eosipova/Brood_parasites/
ref=molAte
query=ageTri
query=agePho

ref_fa=$BP/$ref/Genome/$ref.fa
query_fa=$BP/$query/Genome/$query.fa

minimap2 -cx asm20 --cs $ref_fa $query_fa > $query.$ref.minimap2.paf

# DONE in ~2 hours!
```

#### Call variants from alignment
```
cat $query.$ref.minimap2.paf | sort -k6,6 -k8,8n | paftools.js call -f $ref_fa - | sed "s/sample/$query.minimap/" > $query.$ref.from_minimap2.vcf
```

#### Filter outgroup VCFs for Degenotate: remove indels!
```
bcftools view $query.$ref.from_minimap2.vcf -v snps -U -m2 -M2 -e 'F_MISSING > 0.75 | ref="N" | ALT="." | TYPE~"indel"' -O z -o filt.$query.$ref.from_minimap2.vcf.gz
```



### 1.5. Run Degenotate

[Degenotate source code and installation](https://github.com/harvardinformatics/degenotate)

#### merge 3 VCFs
```
BP=/n/holylfs05/LABS/informatics/Lab/project-eosipova/Brood_parasites/
db=molAte
WDIR=$BP/$db/
cd $WDIR

mA_VCF=Degenotate/filt_strict.molAte.passed_raw.vcf.gz
aP_VCF=Alignments/vs_agePho/filt.agePho.molAte.from_minimap2.vcf.gz
aT_VCF=Alignments/vs_ageTri/filt.ageTri.molAte.from_minimap2.vcf.gz
VCF=Degenotate/filt_strict.mA.aP_aT_from_minimap.vcf.gz

bcftools merge --missing-to-ref $mA_VCF $aP_VCF $aT_VCF -O z -o $VCF
tabix $VCF
```

#### run imputed MK test with Degenotate
```
GENOME=$BP/$db/Genome/$db.fa
ANNO=$BP/$db/Annotation/$db.ncbi_anno.gtf
outgroups=agePho.minimap,ageTri.minimap
DEGENOTATE=/n/holylfs05/LABS/informatics/Lab/project-eosipova/Scripts/degenotate/degenotate.py

conda activate degenotate

$DEGENOTATE -a $ANNO -g $GENOME -u $outgroups -v $VCF
```

#### Get the longest transcript
```
filter_annotation_with_list.py -c 1 -a mk.tsv -l <(awk '$5=="yes"{print $1}' transcript-counts.tsv) > longest.mk.tsv
```

#### Add gene labeles
```
ISO=../../Annotation/ncbi_isoforms.tsv
renameToHLscaffolds.py -c 1 -a longest.mk.tsv -d <(awk '{print $2","$0}' $ISO) > gene.longest.mk.tsv
```


the next steps are performed locally

### 1.6. Analize McDonald-Kreitman test results
```
DBS="vidCha vidMac anoImb indInd molAte poeAcu picPub agePho"
```


#### 1.6.1. Get lists of genes under positive selection & genes with DoS<0 (p<0.01)
```
for db in $DBS; \
do \
	echo $db; \
	awk '($10<0.01 && $12>0){print $1}' MK_test_${db}_ncbi/imp.gene.longest.mk.tsv > MK_test_${db}_ncbi/impMKT/pos.genes.p0_01.lst; \
	awk '($10<0.01 && $12<0){print $1}' MK_test_${db}_ncbi/imp.gene.longest.mk.tsv > MK_test_${db}_ncbi/impMKT/neg.genes.p0_01.lst; \
done
```


#### 1.6.2. GO Enrichment analysis of imputed MK test results
NB: rename bird genes based on the association of chicken gene names with human gene names
```
RENAMEDICT=~/Documents/LabDocs/Chicken_resources/galGal6_gene.hg38_gene_symbol.tsv

for db in $DBS\
do \
	echo $db; \
	
	renameToHLscaffolds.py -c 1 -a MK_test_${db}_ncbi/imp.gene.longest.mk.tsv -d <(sed 's/\t/,/' $RENAMEDICT) | grep -v ^reg_ | grep -v ^LOC > MK_test_${db}_ncbi/hg38.imp.gene.longest.mk.tsv; \

	run_goenrich_analysis.R -w $(pwd) -g  MK_test_${db}_ncbi/hg38.imp.gene.longest.mk.tsv -o MK_test_${db}_ncbi/impMKT/ -p 0.01; \
done
```


#### 1.6.3. Remove children GO terms
```
GOOBO=~/Documents/LabDocs/GO_terms_genes/go.obo
#DOS=pos
DOS=neg
ENRICH=$DOS.enrichGO.all_genes_BG.tsv
c=1

for db in $DBS; \
do \
	echo $db; \
	WDIR=MK_test_${db}_ncbi/impMKT/
	f=$WDIR/$ENRICH; \
	golist=$(cut -f$c $f  | tail +2 | tr '\n' ',')

	for g in $(cut -f$c $f | tail +2); \
	do get_go_children.py -f $GOOBO -go $g -l $golist; done | grep "has parents" | awk '{print $1}' > to_exclude_go.lst; \
	
	filter_annotation_with_list.py -b -c $c -a $f -l to_exclude_go.lst > $WDIR/noChildren.$ENRICH; \
done
```


#### 1.6.4. GSEA analysis of imputed MK test results
```
for db in $DBS; \
do \
	echo $db; \
	cut -f1,12 MK_test_${db}_ncbi/imp.gene.longest.mk.tsv | grep -v -w NA | grep -v ^LOC | tail -n +2 > MK_test_${db}_ncbi/impMKT/gene.dos.tsv; \

	run_gsea_analysis.R -w $(pwd) -g  MK_test_${db}_ncbi/impMKT/gene.dos.tsv  -o MK_test_${db}_ncbi/impMKT/gse.tsv -t MK_test_${db}_ncbi/imp.gene.longest.mk.tsv; \
done
```


#### 1.6.5. Test convergence: make table for GO terms occuring at least in 2 clades
```
for db in $DBS; \
do \
	cut -f1 MK_test_${db}_ncbi/impMKT/pos.enrichGO.all_genes_BG.count3.tsv; \
done | g -v ID | s | uniq -c | awk '$1>1{print $2}' > go_convergent_2clades.v2.lst


f=impMKT/pos.enrichGO.all_genes_BG.count3.tsv

for go in $(cat go_convergent_2clades.v2.lst ); \
do \
	description=$(grep -w $go MK_test_*_ncbi/$f | head -1 | cut -f1,2 | cut -d: -f2-); \
	all_pval=''; \
	
	for db in $DBS; \
	do \
		pval=$(grep -w $go MK_test_${db}_ncbi/$f | cut -f5); \
		if [[ -z $pval ]]; \
		then \
			pval=NA; \
		fi; \
		all_pval=$all_pval"\t"$pval; \
	done; \

	echo -e "$description\t$all_pval"; \
done > go_convergent_2clades.v2.tsv
```

### MKT visualization
Further analysis and visualization can be found in this [jupyter-notebook](https://github.com/osipovarev/Brood_parasites_analysis/blob/main/mk_test_analysis.ipynb)


