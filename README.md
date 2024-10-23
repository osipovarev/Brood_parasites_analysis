
## 0. Get the tree from OpenTree
```
source ~/.bash_profile
source ~/.bashrc

python getFromOpenTree.py -j bp_birds_sci_names.json -c out.cites -t out.tre
```


## 1. Analize McDonald-Kreitman test results
```
DBS="vidCha vidMac anoImb indInd molAte poeAcu picPub agePho"
```


### 1.1. GO Enrichment analysis of imputed MK test results
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


### 1.1.1. Remove children GO terms
```
GOOBO=~/Documents/LabDocs/GO_terms_genes/go.obo

ENRICH=pos.enrichGO.all_genes_BG.tsv
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


### 1.2. GSEA analysis of imputed MK test results
```
for db in $DBS; \
do \
	echo $db; \
	cut -f1,17 MK_test_${db}_ncbi/imp.gene.longest.mk.tsv | grep -v -w NA | grep -v ^LOC | tail -n +2 > 
	MK_test_${db}_ncbi/impMKT/gene.dos.tsv; \

	run_gsea_analysis.R -w $(pwd) -g  MK_test_${db}_ncbi/impMKT/gene.dos.tsv  -o MK_test_${db}_ncbi/impMKT/gse.tsv -t MK_test_${db}_ncbi/imp.gene.longest.mk.tsv; \
done
```


### 1.3. Test convergence: make table for GO terms occuring at least in 2 clades
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



#######################

## 2. SweepFinder2: overlap LR reaks with low PI: GO enrichement
```
for db in $DBS; \
do \
	echo $db; \
	goenrich_genelist.R -w $(pwd) -g PopGen_${db}/genes.SF2_peaks_low_PI_depth100.1Mb_domain.lst -o  PopGen_${db}/SF2.go_enrich.1Mb_domain.tsv -u PopGen_${db}/genes.ncbi.lst; \
done
```


### 2.1. Remove children GO terms from the list of 3-way convergent terms
```
f=SF2.go_convergent_3clades.tsv
c=2

golist=$(cut -f$c $f  | tail +2 | tr '\n' ',')

for g in $(cut -f$c $f | tail +2); \
do \
	get_go_children.py -f $GOOBO -go $g -l $golist; \
done | grep "has parents" | awk '{print $1}' > to_exclude_go.lst
```
### use this list for plotting 3-way convergence! (see GO_enrichment_analysis.ipynb)


### optional: get resulting table with no children GO terms
```
filter_annotation_with_list.py -b -c $c -a $f -l to_exclude_go.lst > noChildren.$f
```



## Do window based PI and TajimaD for reagions with high LR in SweepFinder2

### window = 100kb; step = 1000bp; coordinate of the interval = middle of the 10kb interval
```
for db in $DBS; \
do \

 ## copy
 scp eosipova@boslogin01.rc.fas.harvard.edu://n/holylfs05/LABS/informatics/Lab/project-eosipova/Brood_parasites/$db/PopGen/high_SF2.$db.slide* PopGen_${db}/; \
 
 ## rename PI; assign START to middle of an interval
 renameToHLscaffolds.py -c 1 -a PopGen_${db}/high_SF2.$db.slide.windowed.pi -d <(awk '{print$2","$0}' PopGen_${db}/chrom.scaffold.dict.tsv) | awk '{$3=$3+($4-$3)/2; print }' > file; mv file PopGen_${db}/high_SF2.$db.slide.windowed.pi;

 ## rename TajimaD; assign START to middle of an interval
 renameToHLscaffolds.py -c 1 -a PopGen_${db}/high_SF2.$db.slide.TajimaD -d <(awk '{print$2","$0}' PopGen_${db}/chrom.scaffold.dict.tsv) |  awk '{$3=$3+($4-$3)/2; print }' > file; mv file PopGen_${db}/high_SF2.$db.slide.TajimaD; \

done
```





