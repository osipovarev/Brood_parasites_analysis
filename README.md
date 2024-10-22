
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
```
for db in $DBS\
do \
	echo $db; \
	run_goenrich_analysis.R -w $(pwd) -g  MK_test_${db}_ncbi/imp.gene.longest.mk.tsv -o MK_test_${db}_ncbi/impMKT/ -p 0.01; \
	## add step to filter out Hit counts < 3:
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





