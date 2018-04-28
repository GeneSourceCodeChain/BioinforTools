#!/bin/bash
#
# 2018-4-28
# This script is for download the hg38 version of the human genome template.

for i in $(seq 1 22) X Y M;
do echo $i;
wget ftp://hgdownload.cse.ucsc.edu/goldenPath/hg38/chromosomes/chr${i}.fa.gz;
done
gunzip *.gz
for i in $(seq 1 22) X Y M;
do cat chr${i}.fa >> hg38.fa;
done
rm -fr chr*.fa
