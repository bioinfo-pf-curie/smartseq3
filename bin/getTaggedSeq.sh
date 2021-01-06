#!/bin/bash

fastq_in1=$1
fastq_in2=$2
fastq_out1=$3
txt_out2=$4
fastq_out3=$5

# get tagged sequences
seqkit grep --by-seq --pattern "ATTGCGCAATG" $fastq_in1 > $fastq_out1
echo $fastq_out1
head $fastq_out1
# extract ids
seqkit seq -n -i $fastq_out1 > $txt_out2
echo $txt_out2
head $txt_out2
# get paired reads
seqkit grep -f $txt_out2 $fastq_in2 -o $fastq_out3
echo $fastq_out3
head $fastq_out3
