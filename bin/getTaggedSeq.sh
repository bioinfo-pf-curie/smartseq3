#!/bin/bash

fastq1=$1
fastq2=$3
fastq_out1=$5
txt_out2=$7
fastq_out3=$9


# get tagged sequences
seqkit grep --by-seq --pattern "ATTGCGCAATG" $fastq1 -o $fastq_out1
# extract ids
seqkit seq -n $fastq_out1 > $txt_out2
# get paired reads
seqkit grep --pattern-file $txt_out2 $fastq2 -o $fastq_out3