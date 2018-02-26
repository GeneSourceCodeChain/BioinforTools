#!/bin/bash

workdir=/mnt/harddrive/light/shell/shell_5
for i in aaa bbb ccc ddd eee
do
    mkdir $workdir/${i}_diff
    cp $workdir/diff.sh $workdir/${i}_diff/
    sed -r -i -e "s/$i/aaa/" -e "s/(sample=).*/\1$i/" $workdir/${i}_diff/diff.sh
    /bin/bash $workdir/${i}_diff/diff.sh
    rm -rf $workdir/${i}_diff/diff.sh
done
