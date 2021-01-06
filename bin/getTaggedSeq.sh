#!/bin/bash


# get tagged sequences
seqkit grep --by-seq --pattern "ATTGCGCAATG" $2 > $4
echo $4
head $4
# extract ids
seqkit seq -n -i $4 > $5
echo $5
head $5
# get paired reads
seqkit grep -f $5 $3 -o $6
echo $6
head $6
