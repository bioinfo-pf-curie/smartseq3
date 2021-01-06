#!/bin/bash

echo $2
head $2

echo $3
head $3

echo $4
head $4

echo $5
head $5

echo $6
head $6

# get tagged sequences
seqkit grep --by-seq --pattern "ATTGCGCAATG" $2 > $4
# extract ids
seqkit seq -n -i $4 > $5
# get paired reads
seqkit grep -f $5 $3 -o $6
