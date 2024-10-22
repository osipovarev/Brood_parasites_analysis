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

