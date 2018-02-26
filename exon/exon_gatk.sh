#!/bin/bash
workdir=/the/path/to/dir/S072_12_08
fq1=/the/path/to/dir/sample/S072_S4_L008_R1_001.fastq.gz
fq2=/the/path/to/dir/sample/S072_S4_L008_R2_001.fastq.gz
sample=S072

GENOME=/the/path/to/annotation/hg19/hg19/hg19.fa
INDEX=/the/path/to/annotation/hg19/hg19/hg19.fa
GATK=/the/path/to/dir/biosoft/GATK/GenomeAnalysisTK.jar
PICARD=/the/path/to/dir/biosoft/picardtools/2.9.2/picard.jar
picarddir=/the/path/to/dir/biosoft/picardtools/picard-tools-1.119
DBSNP=/the/path/to/annotation/hg19/bundle/dbsnp_138.hg19.vcf.gz
SNP=/the/path/to/annotation/hg19/bundle/1000G_phase1.snps.high_confidence.hg19.sites.vcf.gz
OM_SITE=/the/path/to/annotation/hg19/bundle/1000G_omni2.5.hg19.sites.vcf
Mills_indels=/the/path/to/annotation/hg19/bundle/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf
KG_indels=/the/path/to/annotation/hg19/bundle/1000G_phase1.indels.hg19.sites.vcf
TMPDIR=/the/path/to/dir/tmp/software
threads=`grep 'processor' /proc/cpuinfo | sort -u | wc -l`
bedfile=/the/path/to/annotation/hg19/hg19_bed/hg19_exons.bed
#######################################################################
# 1 Alignment
echo "start alignment."
echo `date`
bwa mem -t $threads -M -R "@RG\tID:$sample\tSM:$sample\tLB:WES\tPL:Illumina" $INDEX $fq1 $fq2 1>$sample.sam 2>/dev/null

# 2 Sort and Index
echo "start sort and index."
echo `date`
java -Djava.io.tmpdir=$TMPDIR -Xms2g -Xmx10g -jar $PICARD SortSam SORT_ORDER=coordinate INPUT=$sample.sam OUTPUT=$sample.bam
samtools index $sample.bam
rm $sample.sam

# 3 Basic Statistics
echo "start flagstat and idxstats."
echo `date`
samtools flagstat $sample.bam > ${sample}.alignment.flagstat
samtools idxstats $sample.bam > ${sample}.alignment.stat

# 4 multiple filtering for bam files
echo "start marker for bam files."
echo `date`
ulimit -c unlimited \
&& java -Djava.io.tmpdir=$TMPDIR -Xms2g -Xmx10g -jar $picarddir/MarkDuplicates.jar \
INPUT=$sample.bam OUTPUT=${sample}_marker.bam METRICS_FILE=$sample.metrics
rm $sample.bam
echo "start fixed."
echo `date`
java -Djava.io.tmpdir=$TMPDIR -Xms2g -Xmx10g -jar $picarddir/FixMateInformation.jar \
INPUT=${sample}_marker.bam OUTPUT=${sample}_marked_fixed.bam SO=coordinate
samtools index ${sample}_marked_fixed.bam
rm ${sample}_marker.bam

# 5 gatk process bam files
echo "start split."
echo `date`
java -Djava.io.tmpdir=$TMPDIR -Xms2g -Xmx10g -jar $GATK -T SplitNCigarReads \
 -R $GENOME -I ${sample}_marked_fixed.bam -o ${sample}_marked_fixed_split.bam \
 -rf ReassignOneMappingQuality -RMQF 255 -RMQT 60 -U ALLOW_N_CIGAR_READS
rm ${sample}_marked_fixed.bam
echo "start Realigner."
echo `date`
java -Djava.io.tmpdir=$TMPDIR -Xms2g -Xmx10g -jar $GATK -T RealignerTargetCreator \
 -I ${sample}_marked_fixed_split.bam -R $GENOME -o ${sample}_target.intervals \
 -known $Mills_indels -known $KG_indels -known $OM_SITE -nt 5

echo "start IndelRealigner."
echo `date`
ulimit -c unlimited \
&& java -Djava.io.tmpdir=$TMPDIR -Xms2g -Xmx10g -jar $GATK -T IndelRealigner \
 -I ${sample}_marked_fixed_split.bam  -R $GENOME -targetIntervals ${sample}_target.intervals \
 -o ${sample}_realigned.bam
 -known $Mills_indels -known $KG_indels
rm ${sample}_marked_fixed_split.bam
echo "start recall."
echo `date`
java -Djava.io.tmpdir=$TMPDIR -Xms2g -Xmx10g -jar $GATK -T BaseRecalibrator \
 -I ${sample}_realigned.bam -R $GENOME -o ${sample}_temp.table -knownSites $SNP
echo "start print reads."
echo `date`
java -Djava.io.tmpdir=$TMPDIR   -Xms2g -Xmx10g -jar $GATK -T PrintReads \
 -R $GENOME -I ${sample}_realigned.bam -o ${sample}_recal.bam -BQSR ${sample}_temp.table
samtools index ${sample}_recal.bam
rm ${sample}_realigned.bam
chmod a=r ${sample}_recal.bam

# 6 gatk call snp/indel
echo "start haplotype."
echo `date`
java -Djava.io.tmpdir=$TMPDIR   -Xms2g -Xmx10g -jar $GATK -T HaplotypeCaller  \
 -L $bedfile \
 -rf  BadCigar \
 -minReadsPerAlignStart 5 \
 -R $GENOME -I ${sample}_recal.bam --dbsnp $DBSNP  \
 -o  ${sample}_raw.snps.indels.vcf
echo "start select variants."
echo `date`
java -Djava.io.tmpdir=$TMPDIR   -Xms2g -Xmx10g -jar $GATK -T SelectVariants  -R $GENOME  \
 -selectType SNP \
 -V ${sample}_raw.snps.indels.vcf -o ${sample}_raw_snps.vcf
java -Djava.io.tmpdir=$TMPDIR   -Xms2g -Xmx10g -jar $GATK -T SelectVariants  -R $GENOME \
 -selectType INDEL  \
 -V ${sample}_raw.snps.indels.vcf   -o ${sample}_raw_indels.vcf

# for SNP
echo "start snp."
echo `date`
java -Djava.io.tmpdir=$TMPDIR   -Xms2g -Xmx10g -jar $GATK -T VariantFiltration -R $GENOME  \
 --filterExpression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || SOR > 4.0 || ReadPosRankSum < -8.0"  \
 --filterName "my_snp_filter" \
 -V ${sample}_raw_snps.vcf  -o ${sample}_filtered_snps.vcf

java -Djava.io.tmpdir=$TMPDIR   -Xms2g -Xmx10g -jar $GATK -T SelectVariants -R $GENOME  \
 --excludeFiltered \
 -V ${sample}_filtered_snps.vcf  -o  ${sample}_filtered_PASS.snps.vcf

java -Djava.io.tmpdir=$TMPDIR   -Xms2g -Xmx10g -jar $GATK -T VariantEval -R $GENOME  \
 -eval ${sample}_filtered_PASS.snps.vcf -o  ${sample}_filtered_PASS.snps.vcf.eval

# for INDEL
echo "start indel."
echo `date`
java -Djava.io.tmpdir=$TMPDIR   -Xms2g -Xmx10g -jar $GATK -T VariantFiltration -R $GENOME  \
 --filterExpression "QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0 || SOR > 10.0"  \
 --filterName "my_indel_filter" \
 -V ${sample}_raw_indels.vcf  -o ${sample}_filtered_indels.vcf

java -Djava.io.tmpdir=$TMPDIR   -Xms2g -Xmx10g -jar $GATK -T SelectVariants -R $GENOME  \
 --excludeFiltered \
 -V ${sample}_filtered_indels.vcf  -o  ${sample}_filtered_PASS.indels.vcf

java -Djava.io.tmpdir=$TMPDIR   -Xms2g -Xmx10g -jar $GATK -T VariantEval -R $GENOME  \
 -eval ${sample}_filtered_PASS.indels.vcf -o  ${sample}_filtered_PASS.indels.vcf.eval

