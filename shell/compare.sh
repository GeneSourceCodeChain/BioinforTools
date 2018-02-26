#!/bin/bash

# The script is to compare the Y chromosome genotypes of multiple individuals, and to identify individuals with identical genotypes.

workdir=/mnt/harddrive/light/shell/shell_6
sample=$workdir/Y_gene.txt
n=`cat $sample | wc -l`
line=`awk -F "\t" '{print NF}' $sample | sed -n 1p`
monthday=`date +%Y%m%d`
dir1=$workdir/result

mkdir -p $dir1/equal
touch $dir1/equal/all.xlsx

for ((i=1;i<$line;i++))
do
    awk -F "\t" "{print $"$i"}" $sample >$workdir/gene_1.txt
    samp1=$workdir/gene_1.txt
    sam1=`sed -n '1p' $samp1`
    mkdir $workdir/$sam1
    echo "===================================================" >>$dir1/equal/$sam1.txt
    echo "$sam1	identical" >>$dir1/equal/$sam1.txt
    for ((j=`expr $i + 1`;j<=$line;j++))
    do
        awk -F "\t" "{print $"$j"}" $sample >$workdir/gene_2.txt
        samp2=$workdir/gene_2.txt
        sam2=`sed -n '1p' $samp2`
        for k in `seq 2 $n`
        do
            a=`sed -n "${k}p" $samp1`
            b=`sed -n "${k}p" $samp2`
            if [[ $a == $b ]];then
                echo "equality" >>$workdir/$sam1/${sam1}_${sam2}.txt
            elif [[ $a == "" || $b == "" ]];then
                echo "null" >>$workdir/$sam1/${sam1}_${sam2}.txt
            elif [[ $a != $b ]];then
                echo "inequality" >>$workdir/$sam1/${sam1}_${sam2}.txt
            else
                echo "error" >>$workdir/$sam1/${sam1}_${sam2}.txt
            fi
        done
        cat $workdir/$sam1/${sam1}_${sam2}.txt | grep "inequality" >$workdir/a.txt
        if [ $? -eq 0 ];then
            echo "$sam2	different" >>$dir1/$sam1.txt
        else
            echo "$sam2	identical" >>$dir1/$sam1.txt
        fi
    done
    cat $dir1/$sam1.txt | grep "identical" >>$dir1/equal/$sam1.txt
    grep -wFf $dir1/equal/$sam1.txt $dir1/equal/all.xlsx |grep "identical"
    if [ $? -eq 0 ];then
        echo aa >/dev/null
    else
        cat $dir1/equal/$sam1.txt >>$dir1/equal/all.xlsx
    fi
    rm -rf $workdir/$sam1
done

rm -rf $workdir/gene_1.txt $workdir/gene_2.txt $workdir/a.txt

dir2=$workdir/result_$monthday
if [[ -d $dir2 ]];then
    for x in `seq 1 5`
    do
        if [[ -d ${dir2}_$x ]];then
            continue
        else
            mv $dir1 ${dir2}_$x
            break
        fi
    done
else
    mv $dir1 $dir2
fi
