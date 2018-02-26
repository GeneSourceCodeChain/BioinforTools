#!/bin/bash

# Use grep to extract the same line from the beginning string in a file to the new file.

workdir=/mnt/harddrive/light/shell/shell_4

for i in {1..5}
do
    grep -w "^$i" $workdir/test.txt >$workdir/$i.txt
done
