#!/bin/bash

#2018-9-4

workdir=/mnt/harddrive/link_ana/minimac4
dir1=$workdir/panels/m3vcf
dir2=/mnt/harddrive/light/WGC083268D_combined/snv/file/WGC_X_common
resdir=$workdir/result/X_0.05

mkdir -p $resdir

for i in {1..22}
do
    minimac4 --refHaps $dir1/${i}.1000g.Phase3.v5.With.Parameter.Estimates.m3vcf.gz \
	    --haps $dir2/WGC_X_common_chr${i}.vcf \
	    --minRatio 0.05 \
	    --prefix $resdir/X_chr$i &
done
