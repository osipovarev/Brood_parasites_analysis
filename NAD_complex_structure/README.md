
## get anoImb nuclear sequences
```
for i in $(ls ND_alignments/ | cut -d_ -f1 | s -u); do get_seq_by_name_fasta.py -n Aimb -f ND_alignments/${i}_aa.fasta | sed 's/Aimb/Aimb_long/'  > anoImb_mt_ND_models/anoImb.$i.aa.fasta; done
```


## get anoImb MT sequences
```
for i in $(ls mt_alignments/ | cut -d_ -f1 | s -u); do get_seq_by_name_fasta.py -n Aimb_genome -f mt_alignments/${i}_35tax_aa.fasta  | sed 's/Aimb_genome/Aimb_long/'  > anoImb_mt_ND_models/anoImb.$i.aa.fasta; done
```


## Get sites unique in cuckoo finches
```
for i in $(ls ND_alignments/ | cut -d_ -f1 | s -u); \
do \
	sites=$(get_variable_sites_msa.py -f ND_alignments/${i}_aa.fasta -l Aimb | cut -d: -f2); echo -e "$i$sites"; \
done | tr ' ' '\t'  > ND_genes_Aimb.variable_sites.tsv


reseq_Aimb=Aimb_genome,2013ERY045A2,CF0724A,2012JUN006,2014JUN21A1,CF0715B,2014PS058A1,2012PS024A1,2013_PS171A1

for i in $(ls mt_alignments/ | cut -d_ -f1 | s -u); \
do \
	sites=$(get_variable_sites_msa.py -f mt_alignments/${i}_35tax_aa.fasta -l $reseq_Aimb | cut -d: -f2); \
	echo -e "$i$sites";  \
done | tr ' ' '\t'  >> ND_genes_Aimb.variable_sites.tsv
```


## Get sites that are variable in general
```
DBS=Apho,Mate,Scan,Tgut,Lstr,Vmac,Vcha,Aimb

for i in $(ls ND_alignments/ | cut -d_ -f1 | s -u); \
do \
	sites=$(get_variable_sites_msa.py -f ND_alignments/${i}_aa.fasta -l $DBS | cut -d: -f2); echo -e "$i$sites"; \
done | tr ' ' '\t' >  ND_genes.variable_sites.tsv


```


## map sites of interest onto models
```
#VARIABLE=ND_genes_Aimb.variable_sites.tsv
VARIABLE=ND_genes.variable_sites.tsv

for g in $(ls anoImb_mt_ND_models/*pdb | cut -d/ -f2 | cut -d. -f2 | s -u | g ^ND); \
do \
	sites=$(grep -w $g $VARIABLE | cut -f2); \
	if [ -z $sites ]; \
	then \
		cp anoImb_mt_ND_models/anoImb.$g.model.pdb anoImb_mt_ND_models/anoImb.$g.model.BEB.pdb; \
	else map_sites_on_structure.py -p anoImb_mt_ND_models/anoImb.$g.model.pdb -s $sites -f anoImb_mt_ND_models/anoImb.$g.aa.fasta > anoImb_mt_ND_models/anoImb.$g.model.BEB.pdb; \
	fi; \
done
```


## color models by delta dN (repalce B-factor with delta dN values)
```
for g in $(ls anoImb_mt_ND_models/*pdb | cut -d/ -f2 | cut -d. -f2 | s -u); do echo $g; b=$(grep -w $g dN_values.ND_genes.tsv | cut -f5);  awk -v var=$b '{if ($1 == "ATOM" || $1 == "HETATM") printf "%-60s%6.2f%s\n", substr($0, 1, 60), var, substr($0, 67); else print $0}' anoImb_mt_ND_models/anoImb.$g.model.pdb > anoImb_mt_ND_models/anoImb.$g.model.ddN.pdb; done
```


## write pymol script
```
for g in $(ls anoImb_mt_ND_models/*pdb | cut -d/ -f2 | cut -d. -f2 | s -u | g ^ND); do echo "load anoImb_mt_ND_models/anoImb.$g.model.ddN.pdb"; echo "align anoImb.$g.model.ddN, 5lnk "; done > load_and_align_all_NDs.pml
```




## PyMol
```
select beb_nd1, resi 46+79+143+163+168+223+254+258+262+314+315 and anoImb.ND1.model.BEB
select beb_nd3, resi 10+20+62+84+87+94 and anoImb.ND3.model
select bebe_ndufs2, resi 29 and anoImb.NDUFS2.model 

bg_color white
spectrum b, blue_white_red, minimum=-5, maximum=5
hide cartoon, !(ss H+S)

select beb, b>0 and name ca
```

## Run jupyter-notebook => extract shortest distances to each nuclear-encoded residue

## Reverse residue numbering to the original order
```
renameToHLscaffolds.py -c 4 -a <(renameToHLscaffolds.py -c 1 -a distances_complex_I_anoImb.tsv -d <(awk '{print $1","$0}' renum_dict.complex_I_anoImb.tsv) ) -d <(awk '{print $1","$0}' renum_dict.complex_I_anoImb.tsv) > numbered_original.distances_complex_I_anoImb.tsv
```


