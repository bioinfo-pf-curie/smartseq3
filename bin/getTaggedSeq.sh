#!/bin/bash

fastq1=$1
fastq2=$2
fastq_out1="$3"
txt_out2="$4"
fastq_out3="$5"

# get tagged sequences
seqkit grep --by-seq --pattern "TGCGCAATG" $fastq1 -o "$fastq_out1"
# extract ids
seqkit seq -i -n "$fastq_out1" -o "$txt_out2" 
# get paired reads
seqkit grep -f $txt_out2 $fastq2 -o "$fastq_out3"

# test in /data/users/lhadjabe/smartSeq3/smartSeq3_V590/smartseq3/work/55/1d5d0ac19426e866cca14e7ac9c375
# getTaggedSeq.sh V590T17.R1.fastq.gz V590T17.R2.fastq.gz V590T17_tagged_inR1.R1.fastq V590T17_taggedReadIDs_inR1.txt V590T17_tagged_inR1.R2.fastq
# getTaggedSeq.sh V590T17_rest.R2.fastq V590T17.R1.fastq.gz V590T17_tagged_inR2.R2.fastq V590T17_taggedReadIDs_inR2.txt V590T17_tagged_inR2.R1.fastq


# initial script
#seqkit grep --by-seq --pattern "TGCGCAATG" ${reads[0]} -o ${prefix}_tagged_inR1.R1.fastq
# exctract ids
#seqkit seq -n -i ${prefix}_tagged_inR1.R1.fastq -o ${prefix}_taggedReadIDs_inR1.txt
# create R2
#seqkit grep -f ${prefix}_taggedReadIDs_inR1.txt ${reads[1]} -o ${prefix}_tagged_inR1.R2.fastq

#seqkit grep -v -f ${prefix}_taggedReadIDs_inR1.txt ${reads[1]} -o ${prefix}_rest.R2.fastq

# Get tagged sequences in R2 of lefting reads == umi sequences
#seqkit grep --by-seq --pattern "ATTGCGCAATG" ${prefix}_rest.R2.fastq -o ${prefix}_tagged_inR2.R2.fastq 
# exctract ids
#seqkit seq -n -i ${prefix}_tagged_inR2.R2.fastq -o ${prefix}_taggedReadIDs_inR2.txt
# create R1
#seqkit grep -f ${prefix}_taggedReadIDs_inR2.txt ${reads[0]} -o ${prefix}_tagged_inR2.R1.fastq
