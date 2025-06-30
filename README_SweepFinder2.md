

#########################################
#                               		#
##    2. Selective Seweeps analysis     ##
#                               		#
#########################################


db=molAte
BP=/n/holylfs05/LABS/informatics/Lab/project-eosipova/Brood_parasites/


#####################################################

### 2.1. Polarize VCF
```
db=molAte
outdb=agePho
BP=/n/holylfs05/LABS/informatics/Lab/project-eosipova/Brood_parasites/

GENOME=$BP/$db/Genome/$db.fa
OUTGROUP=$BP/$db/Alignments/vs_${outdb}/filt.$outdb.$db.from_minimap2.vcf.gz
Q_VCF=$BP/$db/snpArcher_${db}/passed.${db}_full_raw.vcf.gz
VCF_1OUT=passed_${db}.filt_${outdb}.vcf
```

#### 0) make ancestral fasta from outgroup vcf
```
conda activate broodp
bcftools consensus $OUTGROUP -f $GENOME -o filt.$outdb.$db.from_minimap2.fa
```

#### 1) make VCF with ONE outgroup;
NB: make an unzipped VCF! (the next script doesn't handle gzipped vcf well)
```
bcftools merge --missing-to-ref $Q_VCF $OUTGROUP -O v -o $VCF_1OUT
```

#### 2) infer ancestral allele
```
conda activate vcfdo
vcfdo polarize -i $VCF_1OUT -f filt.$outdb.$db.from_minimap2.fa > polarized.$VCF_1OUT
```

#### split large chroms into chunks (2Mb)
```
cd $BP/$db/Genome/

s=2000000
awk '{print $1"\t1\t"$2}' chrom.sizes > chrom.sizes.bed
awk -v size=$s 'BEGIN{OFS="\t"}{split($3,a,";"); for(i=$2; i<=$3; i+=size) {if (i+size-1 <= $3) print $1, i, (i+size-1); else print $1, i, $3}}' chrom.sizes.bed > split_2Mb.chrom.bed
```

#### prepare input for SweepFinder
```
VCF=$BP/$db/PopGen/polarized.passed_poeAcu.filt_vidMac.vcf.gz
CHROMS=$BP/$db/Genome/split_2Mb.chrom.bed

cd $BP/$db/PopGen/
md split_by_chr_SF2
cd split_by_chr_SF2

for C in $(awk '{print $1":"$2"-"$3}' $CHROMS ); \
do \
  md $C; \
  cat template_sbatch.sh  <(echo -e "tabix $VCF $C --print-header > $C/$C.vcf;") <(echo -e "vcftools --counts2 --derived --vcf $C/$C.vcf --stdout > $C/$C.SF2.input;") > $C/sbatch.$C.vcftools.sh;  \
  echo -e "sbatch $C/sbatch.$C.vcftools.sh";\
done > jobs.vcftools

bash jobs.vcftools
```


#### remove polarized sites with 0 frequency; add header
```
for C in $(awk '{print $1":"$2"-"$3}' $CHROMS ); \
do \
  echo $C; \
  cat $C/$C.SF2.input |awk 'NR<=1 {next} {print $2"\t"$6"\t"$4"\t0"}' | awk '{if ($2 != 0 || $4 != 0) print}' | sed '1i position\tx\tn\tfolded' > $C/header.$C.SF2.input; \
done
```

#### prepare jobs SF2 no rec map
```
for C in $(awk '{print $1":"$2"-"$3}' $CHROMS ); \
do \
  cat template_sbatch.sh <(echo -e "SweepFinder2 -f $C/header.$C.SF2.input $C/$C.SF2.spect;\n SweepFinder2 -lg 1000 $C/header.$C.SF2.input $C/$C.SF2.spect $C/norecmap.header.$C.SF2.out") > $C/sbatch.$C.sweepfinder.norecmap.sh; \
  echo -e "sbatch $C/sbatch.$C.sweepfinder.norecmap.sh"; \
done > jobs.sweepfinder.norecmap

bash jobs.sweepfinder.norecmap
```


#### prepare rec map files
```
RECMAP=$BP/Recombination_maps/all_chrom.recombination_map.$db.bed.gz

for C in $(awk '{print $1":"$2"-"$3}' $CHROMS); \
do \
  chr=$(echo $C | cut -d: -f1); \
  start=$(echo $C | cut -d: -f2 | cut -d- -f1); \
  end=$(echo $C | cut -d: -f2 | cut -d- -f2); \
  cat template_sbatch.sh <(echo -e "fill_in_coord_gaps.py -f <(zcat $RECMAP | awk -v v1=$chr -v v2=$start -v v3=$end '(\$1==v1 && \$2>=v2 && \$3<=v3){print}' | cut -f2,4 | sort -n -k2 | sed \"1i \\$start\\\t0\") > $C/$C.recombination_map.txt")  > $C/sbatch.$C.recomb.sh; \
  echo -e "sbatch $C/sbatch.$C.recomb.sh"; \
done > jobs.recomb_map

bash jobs.recomb_map
```


#### prepare jobs SF2 with rec map
```
for C in $(awk '{print $1":"$2"-"$3}' $CHROMS ); \
do \
  cat template_sbatch.sh <(echo -e "SweepFinder2 -lrg 1000 $C/header.$C.SF2.input $C/$C.SF2.spect $C/$C.recombination_map.txt $C/header.$C.SF2.out") > $C/sbatch.$C.sweepfinder.sh; \
  echo -e "sbatch $C/sbatch.$C.sweepfinder.sh"; \
done > jobs.sweepfinder

bash jobs.sweepfinder
```


### 2. SweepFinder2: overlap LR reaks with low PI: GO enrichement
NB: rename bird genes based on the association of chicken gene names with human gene names
```
RENAMEDICT=~/Documents/LabDocs/Chicken_resources/galGal6_gene.hg38_gene_symbol.tsv

for db in $DBS;\
do \
	echo $db; \
	
	renameToHLscaffolds.py -c 1 -a PopGen_${db}/genes.SF2_peaks_low_PI_depth100.1Mb_domain.lst -d <(sed 's/\t/,/' $RENAMEDICT) | grep -v ^reg_ | grep -v ^LOC > PopGen_${db}/hg38.genes.SF2_peaks_low_PI_depth100.1Mb_domain.lst; \

	goenrich_genelist.R -w $(pwd) -g PopGen_${db}/hg38.genes.SF2_peaks_low_PI_depth100.1Mb_domain.lst -o  PopGen_${db}/SF2.go_enrich.1Mb_domain.tsv -u PopGen_${db}/genes.ncbi.lst; \

	cat PopGen_${db}/SF2.go_enrich.1Mb_domain.tsv | awk -F"\t" '$9>2 {print}' > PopGen_${db}/count3.SF2.go_enrich.1Mb_domain.tsv; \
done
```


### 2.1. Remove children GO terms from the list of 3-way convergent terms
```
GOOBO=~/Documents/LabDocs/GO_terms_genes/go.obo
f=SF2.go_convergent_3clades.tsv
c=2

golist=$(cut -f$c $f  | tail +2 | tr '\n' ',')

for g in $(cut -f$c $f | tail +2); \
do \
	get_go_children.py -f $GOOBO -go $g -l $golist; \
done | grep "has parents" | awk '{print $1}' > to_exclude_go.lst

filter_annotation_with_list.py -b -c $c -a $f -l to_exclude_go.lst > noChildren.$f
```
### use this list for plotting 3-way convergence! (see GO_enrichment_analysis.ipynb)



### 2.2. Remove children and GO terms with hit count <3 in each enrichment analsyis
```
for db in $DBS; \
do \
	echo $db; \
	c=1; \
	f=PopGen_${db}/count3.SF2.go_enrich.1Mb_domain.tsv; \
	golist=$(cut -f$c $f  | tail +2 | tr '\n' ','); \

	for g in $(cut -f$c $f | tail +2); 
	do \
		get_go_children.py -f $GOOBO -go $g -l $golist; \
	done | grep "has parents" | awk '{print $1}' > to_exclude_go.lst; \
	filter_annotation_with_list.py -b -c $c -a $f -l to_exclude_go.lst > noChildren.$db.count3.SF2.go_enrich.1Mb_domain.tsv; \
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





