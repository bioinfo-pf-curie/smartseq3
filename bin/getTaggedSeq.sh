#!/bin/bash

fastq1=$1
fastq2=$2
out1=$3
out2=$4
out3=$5

# get tagged sequences
seqkit grep --by-seq --pattern "ATTGCGCAATG" $fastq1 > $out1
# exctract ids
seqkit seq -n -i $out1 > $out2
# create R2
seqkit grep -f $out2 $fastq2 -o $out3