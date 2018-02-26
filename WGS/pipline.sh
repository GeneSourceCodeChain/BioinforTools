#!/bin/bash

##########################################################################################
# README
# Please make sure the fastq file is *_R1.fastq.gz and *_R2.fastq.gz or you can modify the shell
# The file name: *_R1.fastq.gz *_R2.fastq.gz is the sample name
# The result file is vcf file
# USEAGE
# You should install the procedure "speedseq" "annovar" 
# Put the fastq file into the /the/path/to/WGC083268D_combined
# bash pipline.sh
##########################################################################################

speedseq=/the/path/to/speedseq
annovar=/the/path/to/annovar
sample=WGC083268D_combined
work_dir=/the/path/to/workdir/$sample
monthday=`date +%Y%m%d`
align_dir=$work_dir/align$monthday
snv_dir=$work_dir/snv$monthday
sv_dir=$work_dir/sv$monthday
ref_dir=/the/path/to/ref
fastq_dir=/the/path/to/fastq/file
threads=`grep 'processor' /proc/cpuinfo | sort -u | wc -l`

# configure the environment
function config() {
    echo "start configuring"
    mkdir -p $workdir
    mkdir $work_dir/align$monthday $work_dir/snv$monthday $work_dir/sv$monthday
}

# 1 align
function align() {
    echo "start analysing"
    echo "start mapping,mapping......"
    $speedseq align -t $threads \
     -o $align_dir/$sample \
     -R "@RG\tID:$sample\tSM:$sample\tPL:ILLUMINA" \
     $ref_dir/human_g1k_v37.fasta \
     $fastq_dir/${sample}_R1.fastq.gz \
     $fastq_dir/${sample}_R2.fastq.gz \
     >$align_dir/align.log
}

# 2 vcf calling && unzip
function call() {
    echo "start call snp! snpping......"
    $speedseq var \
     -t $threads \
     -w /the/path/to/speedseq/annotations/ceph18.b37.include.2014-01-15.bed \
     -o $snv_dir/$sample \
     $ref_dir/human_g1k_v37.fasta \
     $align_dir/${sample}.bam \
     >$snv_dir/var.log && \
gunzip $snv_dir/${sample}.vcf.gz
}

# vcf annotation
function vcf_annota() {
    echo "start vcf annotation,annotation......"
    perl $annovar/table_annovar.pl \
     $snv_dir/${sample}.vcf \
     $annovar/humandb/ \
     -buildver hg19 \
     -out $snv_dir/$sample \
     -remove \
     -protocol refGene,avsnp147,1000g2015aug_eas,clinvar_20160302,dbnsfp30a \
     -operation g,f,f,f,f \
     -nastring .  \
     -vcfinput \
     --thread $threads \
     >$snv_dir/annotation.log
}

# 3 Detect SVs && unzip
function detect() {
    echo "start detection of SVs, detecting......"
    $speedseq sv \
     -o $sv_dir/$sample \
     -B $align_dir/${sample}.bam \
     -S $align_dir/${sample}.splitters.bam \
     -D $align_dir/${sample}.discordants.bam \
     -R $ref_dir/human_g1k_v37.fasta \
     -x /the/path/to/speedseq/annotations/ceph18.b37.lumpy.exclude.2014-01-15.bed \
     >$sv_dir/sv.log  && \
gunzip $sv_dir/${sample}.sv.vcf.gz
}

# sv annotation
function sv_annota() {
    echo "start sv.vcf annotation, annotating......"
    perl $annovar/table_annovar.pl \
     $sv_dir/${sample}.sv.vcf \
     $annovar/humandb/ \
     -buildver hg19 \
     -out $sv_dir/${sample}.sv \
     -remove \
     -protocol refGene,cytoBands,avsnp147,1000g2015aug_eas,clinvar_20160302,dbnsfp30a \
     -operation g,r,f,f,f,f \
     -nastring .  \
     -vcfinput \
     --thread $threads \
     >$sv_dir/annosv.log
    echo "Congratulations! All pipline is OK!\n"
}

config
align
call
vcf_annota
detect
sv_annota

