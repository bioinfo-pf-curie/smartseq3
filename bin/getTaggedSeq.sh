#!/bin/bash

echo $2
head <(gzip -cd $2)

echo $3
head <(gzip -cd $3)

# get tagged sequences
seqkit grep --by-seq --pattern "ATTGCGCAATG" $2 -o $4
# extract ids
seqkit seq -n -i $4 -o $5
# get paired reads
seqkit grep -f $5 $3 -o $6

echo $4
head $4

echo $5
head $5

echo $6
head $6