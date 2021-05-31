#!/bin/bash

splan=$1

## Catch sample names
all_samples=$(awk -F, '{print $1}' $splan)

## Add warning if some samples have been removed 
if [[ -f "workflowSummary/warnings.txt" ]]
then
    removed_samples=$(cat workflowSummary/warnings.txt)
    desc="!! WARNING !! <br> $removed_samples" 
    sed -i "s|{desc}|$desc|g" ../../../../assets/multiqcConfig.yaml
fi

#else
#    sed -i "s|{desc}||g" ../../../../assets/multiqcConfig.yaml

# to add in future to write wich samples to bad
#

## Header
#Table:
echo -e "Sample_id,Sample_name,Number_of_reads,Percent_UMIs,Percent_Aligned,Percent_Assigned,Number_of_genes,Number_of_UMIs" > table_mqc.stats
#Bargraph:
echo -e "Sample_id,Sample_name,Aligned_Assigned,Aligned_NotAssigned,NotAligned_NotAssagned" > final_mqc.stats

for sample in $all_samples
do
    ##id
    sname=$(awk -F, -v sname=$sample '$1==sname{print $2}' $splan)

    # Number of FRAGMENTS (tot_reads = nb_frgmts * 2) :
    tot_reads=`grep "totReads: " totReads/${sample}_totReads.txt | sed -e 's/totReads: //'`
    pUMIs=`grep "percentUMI:" pUMIs/${sample}_pUMIs.txt | cut -d':' -f2`
    aligned=`grep "Uniquely mapped reads number" star/${sample}Log.final.out | awk '{print $NF}'`
    paligned=`grep "Uniquely mapped reads %" star/${sample}Log.final.out | awk '{print $NF}' | sed -e 's/%//'`
    aligned_assigned=`grep "Assigned" FC/${sample}_counts.summary | cut -f2`
    NotAligned=`echo $(( $tot_reads*2 - $aligned ))`
    aligned_NotAssigned=`echo $(( $aligned - $aligned_assigned ))`

    if [ $aligned_assigned != 0 ]; then  
        paligned_assigned=$(echo "scale=2; ($aligned_assigned*100/($tot_reads*2))" | bc -l)
    else
	    paligned_assigned=0
    fi

    nbGenes=$(grep ${sample}\" resume.txt | cut -d, -f3)
    nbUMIs=$(grep ${sample}\" resume.txt | cut -d, -f2)

    echo -e ${sample},${sname},${tot_reads},${pUMIs},${paligned},${paligned_assigned},${nbGenes},${nbUMIs} >> table_mqc.stats
    echo -e ${sample},${sname},${aligned_assigned},${aligned_NotAssigned},${NotAligned} >> final_mqc.stats

done
