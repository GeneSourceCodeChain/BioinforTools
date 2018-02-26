#!/bin/bash

# Merge multiple files into a file in order.

workdir=/mnt/harddrive/light/shell/shell_2

for i in {1..5}
do
    cat $workdir/$i.txt >>$workdir/comb.txt
done
