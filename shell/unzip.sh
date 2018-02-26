#!/bin/bash

# There are multiple compressed files in the directory, and the file names contained in each compressed file are the same and the contents are different.
# At the same time decompress multiple compression files and rename the extracted files.

workdir=/mnt/harddrive/light/shell_1

for i in {1..5}
do
    unzip $workdir/$i.zip -d $workdir/ \
    && mv $workdir/a.txt $workdir/$i.txt 
done
