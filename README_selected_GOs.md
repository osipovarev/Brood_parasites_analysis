# This README describes analysis of selected GO terms


## 1. Select GO terms that are the widest for the the following processes
```
cat selected_GOs.lst
```
GO:0007283	spermatogenesis
GO:0002376	Immune system process
GO:0007399	Nervous system development


## 2. Identify genes associated with these GO terms
in chicken GO terms

```
GO_DB=validated.final.isoformes.galGal6.ncbi_appris_go.Maggie.tsv

for go in $(cut -f1 selected_GOs.lst); \
do \
	grep -w "$go" $GO_DB | cut -f1 | s -u > $go.genes.galGal6.lst; \
done
```


## 3. Calculate length of these genes
```
BED=validated.final.galGal6.ncbi_appris.anno.bed

for go in $(cut -f1 selected_GOs.lst); \
do \
	echo $go; \
	for g in $(cat $go.genes.galGal6.lst); \
	do  \
		grep $g $BED | head -1 | cut -f11 | tr ',' '\n' | sum_stdin.py; \
	done > $go.gene_lenghts.lst; \
done
```


## 4. Plot length distribution
see jupyter-notebook `GO_enrichment_analysis.ipynb`
