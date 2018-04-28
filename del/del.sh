#!/bin/bash
#
# The purpose of this script is to delete the file size 0 in the current directory.

for filename in `ls`
do
    if test -d $filename
    then b=0
    else    
       a=$(ls -l $filename | awk '{ print $5 }')
            if test $a -eq 0
             then rm $filename
             fi
        fi      
done
