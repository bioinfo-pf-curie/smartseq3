#!/bin/bash

splan=$1

## Catch sample names
all_samples=$(awk -F, '{print $1}' $splan)

## Add warning if some samples have been removed 
if [[ -f "workflowSummary/warnings.txt" ]]
then
    removed_samples=$(cat workflowSummary/warnings.txt)
    sed -i "s|{desc}|$removed_samples|g" multiqcConfig.yaml
else
    sed -i "s|{desc}||g" multiqcConfig.yaml
fi

## Header
#Table:
echo -e "Sample_id,Sample_name,Number_of_fragments,Percent_UMIs,Percent_Aligned,Percent_Assigned,Number_of_UMIs,Number_of_genes,Number_of_reads,Number_of_genes" > table_mqc.stats
#Bargraph:
echo -e "Sample_id,Sample_name,Aligned_Assigned,Aligned_NotAssigned,NotAligned_NotAssagned" > final_mqc.stats

for sample in $all_samples
do
    ##id
    sname=$(awk -F, -v sname=$sample '$1==sname{print $2}' $splan)

    # Number of FRAGMENTS (tot_frag = nb_frgmts * 2) :
    tot_frag=`grep "totFrag: " totFrag/${sample}_nbTotFrag.txt | sed -e 's/totFrag: //'`
    pUMIs=`grep "percentUMI:" pUMIs/${sample}_pUMIs.txt | cut -d':' -f2`
    # star = en frag 
    aligned=`grep "Uniquely mapped reads number" star/${sample}Log.final.out | awk '{print $NF}'`
    paligned=`grep "Uniquely mapped reads %" star/${sample}Log.final.out | awk '{print $NF}' | sed -e 's/%//'`
    # featureCount = assigned en frag
    aligned_assigned=`grep "Assigned" FC/${sample}_counts.summary | cut -f2`
    NotAligned=`echo $(( $tot_frag - $aligned ))`
    aligned_NotAssigned=`echo $(( $aligned - $aligned_assigned ))`

    if [ $aligned_assigned != 0 ]; then  
	paligned_assigned=$(echo "${aligned_assigned} ${tot_frag}" | awk ' { printf "%.*f", 2, $1*100/$2 } ')
    else
	paligned_assigned=0
    fi

    nbUMIs=$(grep ${sample}\" UMI_gene_per_cell.txt | cut -d, -f2)
    nbGenesUmi=$(grep ${sample}\" UMI_gene_per_cell.txt | cut -d, -f3)
    
    nbReads=$(grep ${sample}\" read_gene_per_cell.txt | cut -d, -f2)
    nbGenesRead=$(grep ${sample}\" read_gene_per_cell.txt | cut -d, -f3)

    echo -e ${sample},${sname},${tot_frag},${pUMIs},${paligned},${paligned_assigned},${nbUMIs},${nbGenesUmi},${nbReads},${nbGenesRead}  >> table_mqc.stats
    echo -e ${sample},${sname},${aligned_assigned},${aligned_NotAssigned},${NotAligned} >> final_mqc.stats

done
