#!/bin/bash

# Find a file that is different from other files, and write the results to the new file.

a=0
sample=aaa
workdir=/mnt/harddrive/light/shell/shell_5
sampdir=$workdir/${sample}_diff

# mkdir $workdir/${sample}_diff
cp $workdir/$sample.txt $sampdir/${sample}_1.txt
for i in bbb ccc ddd eee
do
    a=$(($a+1))
    grep -vFf $workdir/$i.txt $sampdir/${sample}_$a.txt >>$sampdir/${sample}_`expr $a + 1`.txt
    rm -rf $sampdir/${sample}_$a.txt
done

