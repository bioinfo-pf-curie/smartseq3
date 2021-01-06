#!/bin/bash

echo $1
head <(gzip -cd $1)

echo $2
head <(gzip -cd $2)

# get tagged sequences
seqkit grep --by-seq --pattern "ATTGCGCAATG" $2 > $4
# extract ids
seqkit seq -n -i $4 > $5
# get paired reads
seqkit grep -f $5 $3 -o $6

echo $3
head $3

echo $4
head $4

echo $5
head $5