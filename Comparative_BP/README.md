## Get a preliminary phylogeny from OpenTree
```
 python /Users/osipova/Documents/Scripts//tree_tools/getFromOpenTree.py  -j main.json -c out.cites -t out.tre
```

## Filter files with inactivating mutations for selected projections only
(on delta)
```
for db in $(cat dbs_comp_bp.lst); \
do \
	echo $db;

	## define files to work with
	genefile=genes.F_2.RR_0.5.lst;
	ORTHO=vs_${db}/orthology_classification.tsv;
	TEMP=vs_${db}/temp/orthology_scores.tsv;
	INACT=vs_${db}/inact_mut_data.txt;
	SELECT_INACT=Selected_inact_mut/selected.inact_mut_data.$db.txt;

	## select projections
	for g in $(cat $genefile); \
	do 
		status=$(grep -w $g $ORTHO | cut -f5); \
		projection=$(grep -w $g $ORTHO | cut -f4); \
		if [[ $status == "one2zero" ]]; \
		then \
			projection=$(grep $g $TEMP | sort -n -k3,3 | tail -1 | awk '{print $1"."$2}'); \
		fi; \
		echo $projection | tr ' ' '\n'; \
	done > selected_projections.lst

	## filer inact_mut file for selected projections
	filter_annotation_with_list.py -c 1 -a <(grep ^# $INACT | sed 's/# //' | sed "s/\t/./") -l selected_projections.lst > $SELECT_INACT; \
done
```
