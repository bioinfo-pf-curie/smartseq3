#!/bin/bash


# get tagged sequences
seqkit grep --by-seq --pattern "ATTGCGCAATG" $1 > $3
echo $3
head $3
# extract ids
seqkit seq -n -i $3 > $4
echo $4
head $4
# get paired reads
seqkit grep -f $4 $2 -o $5
echo $5
head $5
