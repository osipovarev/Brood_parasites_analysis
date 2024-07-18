
## Get the tree from OpenTree
source ~/.bash_profile
source ~/.bashrc

python getFromOpenTree.py -j bp_birds_sci_names.json -c out.cites -t out.tre



## Filter and analize MK test results


### GSEA analysis of MK test results
```
WDIR=$(pwd)
for db in indInd vidCha vidMac; do cut -f1,17 MK_test_${db}_ncbi/imp.gene.longest.mk.tsv | grep -v -w NA | grep -v ^LOC | tail -n +2 > MK_test_${db}_ncbi/impMKT/gene.dos.tsv; run_gsea_analysis.R -w $(pwd) -g  MK_test_${db}_ncbi/impMKT/gene.dos.tsv  -o MK_test_${db}_ncbi/impMKT/gse.tsv -t MK_test_${db}_ncbi/imp.gene.longest.mk.tsv; done
```

### GO Enrichment analysis of MK test results
```
for db in indInd vidCha vidMac; \
do \
	echo $db; \
	run_goenrich_analysis.R -w $(pwd) -g  MK_test_${db}_ncbi/imp.gene.longest.mk.tsv -o MK_test_${db}_ncbi/impMKT/ -p 0.01; \
done
```

### GO enrich and GSEA for imputed MKT
```
for db in vidCha vidMac agePho molAte picPub; \
do \
	echo $db; \
	cut -f1,12 MK_test_${db}_ncbi/imp.gene.longest.mk.tsv | grep -v -w NA | grep -v ^LOC | tail -n +2 > MK_test_${db}_ncbi/impMKT/gene.dos.tsv; \
	run_gsea_analysis.R -w $(pwd) -g  MK_test_${db}_ncbi/impMKT/gene.dos.tsv  -o MK_test_${db}_ncbi/impMKT/gse.tsv -t MK_test_${db}_ncbi/imp.gene.longest.mk.tsv; \
	run_goenrich_analysis.R -w $(pwd) -g  MK_test_${db}_ncbi/imp.gene.longest.mk.tsv -o MK_test_${db}_ncbi/impMKT/ -p 0.05; \
done
```


### make table for GO terms occuring at least in 2 clades
```
for db in vidMac vidCha anoImb poeAcu indInd picPub molAte agePho; do cut -f1 MK_test_${db}_ncbi/impMKT/pos.enrichGO.all_genes_BG.count3.tsv; done | g -v ID | s | uniq -c | awk '$1>1{print $2}' > go_convergent_2clades.v2.lst


f=impMKT/pos.enrichGO.all_genes_BG.count3.tsv

for go in $(cat go_convergent_2clades.v2.lst ); do description=$(grep -w $go MK_test_*_ncbi/$f | head -1 | cut -f1,2 | cut -d: -f2-); all_pval=''; for db in vidMac vidCha anoImb poeAcu indInd picPub molAte agePho; do pval=$(grep -w $go MK_test_${db}_ncbi/$f | cut -f5); if [[ -z $pval ]]; then pval=NA; fi; all_pval=$all_pval"\t"$pval; done; echo -e "$description\t$all_pval"; done > go_convergent_2clades.v2.tsv
```




## PCA honeyguide chr by chr
```
for f in $(ls split_by_chr_PCA_indInd/*eigenvec); do renameToHLscaffolds.py -c 2 -a $f -d  <(awk '{print $1","$0}' snpArcher_QC/iInd.samples.dict) > $f.labeled; done
```
### rename scaffolds to chromosoms
```
for f in $(ls split_by_chr_PCA_indInd/*passed.indInd.pca.eigenvec.labeled); do scaffold=$(echo $f | grep -Eo "CM0[0-9]{5}" );  chr=$(grep $scaffold chroms.sizes | cut -f1); echo -e "mv split_by_chr_PCA_indInd/$scaffold.passed.indInd.pca.eigenvec.labeled split_by_chr_PCA_indInd/$chr.passed.indInd.pca.eigenvec.labeled"; done > rename.sh
```




## Do window based PI and TajimaD for reagions with high LR in SweepFinder2

### window = 100kb; step = 1000bp; coordinate of the interval = middle of the 10kb interval
```
for db in indInd vidCha vidMac; \
do \

 ## copy
 scp eosipova@boslogin01.rc.fas.harvard.edu://n/holylfs05/LABS/informatics/Lab/project-eosipova/Brood_parasites/$db/PopGen/high_SF2.$db.slide* PopGen_${db}/; \
 
 ## rename PI; assign START to middle of an interval
 renameToHLscaffolds.py -c 1 -a PopGen_${db}/high_SF2.$db.slide.windowed.pi -d <(awk '{print$2","$0}' PopGen_${db}/chrom.scaffold.dict.tsv) | awk '{$3=$3+($4-$3)/2; print }' > file; mv file PopGen_${db}/high_SF2.$db.slide.windowed.pi;

 ## rename TajimaD; assign START to middle of an interval
 renameToHLscaffolds.py -c 1 -a PopGen_${db}/high_SF2.$db.slide.TajimaD -d <(awk '{print$2","$0}' PopGen_${db}/chrom.scaffold.dict.tsv) |  awk '{$3=$3+($4-$3)/2; print }' > file; mv file PopGen_${db}/high_SF2.$db.slide.TajimaD; \

done
```

###############################

## OXPHOS genes analysis ####

1) make OXPHOS_gene_list.bed
2) sshD. HLtaeGut5/TOGA/
3) 

```
LIST=autosomal_control_gene_list
LIST=OXPHOS_gene_list

for g in $(cut -f4 $LIST.bed); do line=$(grep -w $g $LIST.bed); for t in $(grep -w $g toga.isoforms.tsv | cut -f2); do echo -e "$t\t$line"; done ; done > $LIST.tsv

for db in HLvidMac2 HLvidCha2 HLanoImb2; do for t in $(cut -f1 autosomal_control_gene_list.tsv); do grep -w $t vs_${db}/orthology_classification.tsv |  cut -f4,5 | awk -v var=$db '{print var"\t"$0}' ; done; done > file

for db in vidMac vidCha anoImb indInd molAte agePho; do grep $db all_bp.oxphos.toga_ncbi.tsv | cut -f2- > MK_test_${db}_ncbi/oxphos.toga_ncbi.tsv; done

for db in vidMac vidCha anoImb indInd molAte agePho; do renameToHLscaffolds.py -c 1 -a MK_test_${db}_ncbi/oxphos.toga_ncbi.dos.tsv -d <(cut -f1,5 OXPHOS_gene_list.tsv | sed "s/\.[0-9]*\t/\t/" | awk '{print $1","$2"\t"$1}') > MK_test_${db}_ncbi/gene.oxphos.toga_ncbi.tsv; done

for db in vidMac vidCha anoImb indInd molAte agePho; do awk -v var=$db '{print var"\t"$0 }' MK_test_${db}_ncbi/gene.oxphos.toga_ncbi.tsv; done > all_bp.gene.oxphos.toga_ncbi.tsv

# prepare table to run ABC-MKT
for db in vidMac vidCha anoImb; do for g in $(cut -f4 autosomal_control_gene_list.bed); do grep -w $g MK_test_${db}_ncbi/extended.af.gene.longest.mk.tsv | awk 'BEGIN {FS="\t"; OFS="\t"} {print $1,$5,","$3,$6,","$4,$7,$8}' | sed 's/|/,/g'; done > MK_test_${db}_ncbi/for_abc.gene.autosomal_control.tsv; done
```



#######################

## SweepFinder2 peaks overlapping low PI: GO enrichement
```
for db in vidMac vidCha indInd molAte agePho picPub anoImb; \
do \
	echo $db; \
	goenrich_genelist.R -w $(pwd) -g PopGen_${db}/genes.SF2_peaks_low_PI_depth100.1Mb_domain.lst -o  PopGen_${db}/SF2.go_enrich.1Mb_domain.tsv -u PopGen_${db}/genes.ncbi.lst; \
done
```


