
## PCA honeyguide chr by chr
```
for f in $(ls split_by_chr_PCA_indInd/*eigenvec); \
do \
	renameToHLscaffolds.py -c 2 -a $f -d  <(awk '{print $1","$0}' snpArcher_QC/iInd.samples.dict) > $f.labeled; \
done
```

### rename scaffolds to chromosoms
```
for f in $(ls split_by_chr_PCA_indInd/*passed.indInd.pca.eigenvec.labeled); do scaffold=$(echo $f | grep -Eo "CM0[0-9]{5}" );  chr=$(grep $scaffold chroms.sizes | cut -f1); echo -e "mv split_by_chr_PCA_indInd/$scaffold.passed.indInd.pca.eigenvec.labeled split_by_chr_PCA_indInd/$chr.passed.indInd.pca.eigenvec.labeled"; done > rename.sh
```

