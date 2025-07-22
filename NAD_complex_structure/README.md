# Structural analysis cuckoo finch respiratory complex I

This pipeline describes structural analysis of the cuckoo finch respiratory complex I
for Sorenson et al. manuscript: "*Accelerated sex chromosome degeneration and mitonuclear coevolution in a brood parasitic bird*"

## 1. Get cucko finch sequences

### 1.1. get anoImb nuclear sequences
```
for i in $(ls ND_alignments/ | cut -d_ -f1 | sort -u); \
do \
	get_seq_by_name_fasta.py -n Aimb -f ND_alignments/${i}_aa.fasta | sed 's/Aimb/Aimb_long/'  > anoImb_mt_ND_models/anoImb.$i.aa.fasta; \
done
```
`get_seq_by_name_fasta.py` script can be found [here](https://github.com/osipovarev/fasta_tools/blob/main/get_seq_by_name_fasta.py)


### 1.2. get anoImb MT sequences
```
for i in $(ls mt_alignments/ | cut -d_ -f1 | sort -u); \
do \
	get_seq_by_name_fasta.py -n Aimb_genome -f mt_alignments/${i}_35tax_aa.fasta  | sed 's/Aimb_genome/Aimb_long/' > anoImb_mt_ND_models/anoImb.$i.aa.fasta; \
done
```



## 2. Get sites unique in cuckoo finches
```
for i in $(ls ND_alignments/ | cut -d_ -f1 | sort -u); \
do \
	sites=$(get_variable_sites_msa.py -f ND_alignments/${i}_aa.fasta -l Aimb | cut -d: -f2); echo -e "$i$sites"; \
done | tr ' ' '\t'  > ND_genes_Aimb.variable_sites.tsv


reseq_Aimb=Aimb_genome,2013ERY045A2,CF0724A,2012JUN006,2014JUN21A1,CF0715B,2014PS058A1,2012PS024A1,2013_PS171A1

for i in $(ls mt_alignments/ | cut -d_ -f1 | sort -u); \
do \
	sites=$(get_variable_sites_msa.py -f mt_alignments/${i}_35tax_aa.fasta -l $reseq_Aimb | cut -d: -f2); \
	echo -e "$i$sites";  \
done | tr ' ' '\t'  >> ND_genes_Aimb.variable_sites.tsv
```
`get_variable_sites_msa.py` script can be found [here](https://github.com/osipovarev/msa_tools/blob/main/get_variable_sites_msa.py)



## 3. Get sites that are variable in general
```
DBS=Apho,Mate,Scan,Tgut,Lstr,Vmac,Vcha,Aimb

for i in $(ls ND_alignments/ | cut -d_ -f1 | sort -u); \
do \
	sites=$(get_variable_sites_msa.py -f ND_alignments/${i}_aa.fasta -l $DBS | cut -d: -f2); echo -e "$i$sites"; \
done | tr ' ' '\t' >  ND_genes.variable_sites.tsv


```


## 4. Model cuckoo finch respiratory complex I subunits from sequences
used https://swissmodel.expasy.org/interactive
and submited sequences of all subunits one by one.
Downloaded resulting pdb models.

## 5. Optional (for PyMol only): map sites of interest onto models
```
#VARIABLE=ND_genes_Aimb.variable_sites.tsv
VARIABLE=ND_genes.variable_sites.tsv

for g in $(ls anoImb_mt_ND_models/*pdb | cut -d/ -f2 | cut -d. -f2 | sort -u | g ^ND); \
do \
	sites=$(grep -w $g $VARIABLE | cut -f2); \
	if [ -z $sites ]; \
	then \
		cp anoImb_mt_ND_models/anoImb.$g.model.pdb anoImb_mt_ND_models/anoImb.$g.model.BEB.pdb; \
	else map_sites_on_structure.py -p anoImb_mt_ND_models/anoImb.$g.model.pdb -s $sites -f anoImb_mt_ND_models/anoImb.$g.aa.fasta > anoImb_mt_ND_models/anoImb.$g.model.BEB.pdb; \
	fi; \
done
```


## 6. Optional (for PyMol only): color models by delta dN (repalce B-factor with delta dN values)
```
for g in $(ls anoImb_mt_ND_models/*pdb | cut -d/ -f2 | cut -d. -f2 | sort -u); \
do \
	echo $g; \
	b=$(grep -w $g dN_values.ND_genes.tsv | cut -f5);  \
	awk -v var=$b '{if ($1 == "ATOM" || $1 == "HETATM") printf "%-60s%6.2f%s\n", substr($0, 1, 60), var, substr($0, 67); else print $0}' anoImb_mt_ND_models/anoImb.$g.model.pdb > anoImb_mt_ND_models/anoImb.$g.model.ddN.pdb; \
done
```
models for individual subunits are available on Dryad



## 7. Load and align models with PyMol

### 7.1. write pymol script
```
for g in $(ls anoImb_mt_ND_models/*pdb | cut -d/ -f2 | cut -d. -f2 | sort -u | g ^ND); \
do \
	echo "load anoImb_mt_ND_models/anoImb.$g.model.ddN.pdb"; echo "align anoImb.$g.model.ddN, 5lnk "; \
done > load_and_align_all_ddN_NDs.pml
```

### 7.2. Run the script in PyMol
=> saved aligned structures as renum.complex_I_anoImb_M_N_ddN.pdb (available on Dryad).


## 8. Compute distance matrix between all mt-encoded residues to nuclear-encoded resudies
Run [jupyter-notebook](https://github.com/osipovarev/Brood_parasites_analysis/blob/main/NAD_complex_analysis.ipynb) 
=> extract shortest distances to each nuclear-encoded residue



## 9. Reverse residue numbering to the original order
```
renameToHLscaffolds.py -c 4 -a <(renameToHLscaffolds.py -c 1 -a distances_complex_I_anoImb.tsv -d <(awk '{print $1","$0}' renum_dict.complex_I_anoImb.tsv) ) -d <(awk '{print $1","$0}' renum_dict.complex_I_anoImb.tsv) > numbered_original.distances_complex_I_anoImb.tsv
```



## 10. Visualize structure in ChimeraX
used [custom_color_resi.pml](https://github.com/osipovarev/Brood_parasites_analysis/blob/main/NAD_complex_structure/custom_color_resi.pml) script to color the resulting structures in ChimeraX.

Saved the ChimeraX session with the structure as NAD_CF_all_CF1_surface_mt.cxs (available on Dryad)

