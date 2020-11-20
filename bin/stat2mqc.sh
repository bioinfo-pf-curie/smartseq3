#!/bin/bash

splan=$1
filtreUmi1=$3
filtreUmi2=$5
filtreGene1=$7
filtreGene2=$9

## Catch sample names
all_samples=$(awk -F, '{print $1}' $splan)


## Header
#Table:
echo -e "Sample_id, Sample_name, Reads tot., UMIs(%), Aligned(%), Assigned(%), Cells tot." > table_mqc.stats
#Bargraph:
#echo -e "Sample_id, Sample_name, Barcoded_Aligned_Assigned , Barcoded_Aligned_NotAssigned " > final_mqc.stats

for sample in $all_samples
do
    ##id
    sname=$(awk -F, -v sname=$sample '$1==sname{print $2}' $splan)

    # Number of reads:
    tot_reads=`grep "Number of input reads" star/${sample}Log.final.out | cut -d'|' -f2`
    pUMIs=`grep "percentUMI:" mergeReads/${sample}_pUMIs.txt | cut -d':' -f2`
    aligned=` grep "Uniquely mapped reads number" star/${sample}Log.final.out | cut -d'|' -f2 `
    aligned_assigned=`grep "Assigned" FC/${sample}_counts.summary | cut -f2`

    # Number of assigned across barcoded and aligned
    NotAligned=`echo $(( $tot_reads - $aligned ))`
    aligned_NotAssigned=`echo $(( $aligned - $aligned_assigned ))`
    NotAligned_NotAssigned=`echo $(( $tot_reads - $aligned_assigned ))`

    if [ $aligned != 0 ]; then  
        pbarcoded_aligned=$(echo "scale=2; ($aligned*100/$tot_reads)" | bc -l)
    else
	    pbarcoded_aligned= 0
    fi

    if [ $aligned_assigned != 0 ]; then  
        paligned_assigned=$(echo "scale=2; ($aligned_assigned*100/$tot_reads)" | bc -l)
    else
	    pbarcoded_aligned_assigned=0
    fi

    ##### Filtre Count
    tot_cells=`wc -l $all_samples`
    #filtre1=`sed -n 2p countsFiltre/${sample}_countsFiltre.log`
    #filtre2=`sed -n 3p countsFiltre/${sample}_countsFiltre.log`
    #filtreG1=`sed -n 4p countsFiltre/${sample}_countsFiltre.log`
    #filtreG2=`sed -n 5p countsFiltre/${sample}_countsFiltre.log`

    echo -e ${sample},${sname},${tot_reads}, ${pUMIs}, ${paligned}, ${paligned_assigned}, ${tot_cells}>> table_mqc.stats
    #echo -e ${sample},${sname},${barcoded_aligned_assigned}, ${barcoded_aligned_NotAssigned},${barcoded_NOTaligned},${NOTbarcoded_aligned},${NOTbarcoded_NOTaligned} >> final_mqc.stats

done




    