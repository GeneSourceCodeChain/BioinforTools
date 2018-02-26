#!/bin/bash

# Merge multiple files with paste.
# Append the next file to the back of the previous file rather than below.

workdir=/mnt/harddrive/light/shell/shell_3

touch $workdir/0.xlsx

for i in {1..5}
do
    j=`expr $i - 1`
    paste $workdir/$j.xlsx $workdir/$i.txt >$workdir/$i.xlsx
    rm -rf $workdir/$j.xlsx
done

sed -i 's/	//' $workdir/5.xlsx
