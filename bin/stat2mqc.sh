#!/bin/bash

splan=$1

## Catch sample names
all_samples=$(awk -F, '{print $1}' $splan)

## Header
#Table:
echo -e "Sample_id, Sample_name, #Fragments, %UMIs, %Aligned, %Assigned, #Genes, #UMIs" > table_mqc.stats
#Bargraph:
echo -e "Sample_id, Sample_name, Aligned_Assigned, Aligned_NotAssigned, NotAligned_NotAssagned" > final_mqc.stats

for sample in $all_samples
do
    ##id
    sname=$(awk -F, -v sname=$sample '$1==sname{print $2}' $splan)

    # Number of reads:
    tot_reads=`grep "Input Reads: " umiExtract/${sample}_umiExtract.log | cut -d":" -f4`
    pUMIs=`grep "percentUMI:" mergeReads/${sample}_pUMIs.txt | cut -d':' -f2`
    aligned=` grep "Uniquely mapped reads number" star/${sample}Log.final.out | cut -d'|' -f2 `
    paligned=` grep "Uniquely mapped reads %" star/${sample}Log.final.out | cut -d'|' -f2 | cut -d% -f1`
    aligned_assigned=`grep "Assigned" FC/${sample}_counts.summary | cut -f2`
    NotAligned=`echo $(( $tot_reads - $aligned ))`
    aligned_NotAssigned=`echo $(( $aligned - $aligned_assigned ))`

    if [ $aligned_assigned != 0 ]; then  
        paligned_assigned=$(echo "scale=2; ($aligned_assigned*100/$tot_reads)" | bc -l)
    else
	    paligned_assigned=0
    fi

    nbGenes=$(grep ${sample} resume.txt | cut -d, -f3)
    nbUMIs=$(grep ${sample} resume.txt | cut -d, -f2)

    echo -e ${sample},${sname},${tot_reads}, ${pUMIs}, ${paligned}, ${paligned_assigned},${nbGenes}, ${nbUMIs} >> table_mqc.stats
    echo -e ${sample},${sname},${aligned_assigned}, ${aligned_NotAssigned},${NotAligned} >> final_mqc.stats

done




    