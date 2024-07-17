#!/usr/bin/env bash

#-------------------------------------------------------------------
# Version:  23.09.30
# Author:   Tan yuanqing
# Function: polyA tail XR-seq analysis
# Usage:    bash XR-polyA-seq_analysis_for_human_hg38.sh Read1 cores
#-------------------------------------------------------------------

# 设置所需变量接受参数
name=$(basename $1 ".fastq.gz")
ref=~/Tanyuanqing/Genome/hg38/bwa_index/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna
chrNamelength=~/Tanyuanqing/Genome/hg38/bwa_index/hg38chrNamelength.txt
cores=$2

# 去除polyA tail 和 adapt,并过滤低测序质量的reads
cutadapt -a "A{10}" -q 30 -j $cores -o ${name}_filterpolyA.fastq $1  
cutadapt -a GATCGGAAGA -q 30 -j $cores -o ${name}_filter.fastq ${name}_filterpolyA.fastq

# short sequence mapping by bwa backtracking algorithm 
bwa aln -t $cores $ref ${name}_filter.fastq > ${name}_filter.sai
bwa samse $ref ${name}_filter.sai ${name}_filter.fastq > ${name}_filter.sam

# sam into bam, sort, remove duplicated
samtools view -b -S -q 20 -@ $cores ${name}_filter.sam > ${name}_filter.bam
sambamba sort -t $cores -o ${name}_sort_filter.bam ${name}_filter.bam
sambamba markdup -r -t $cores ${name}_sort_filter.bam ${name}_markdupremove_sort_filter.bam

# bam into bed, sort
bedtools bamtobed -i ${name}_markdupremove_sort_filter.bam > ${name}_markdupremove_sort_filter.bed
sort -u -k1,1 -k2,2n -k3,3n -k6,6 ${name}_markdupremove_sort_filter.bed > ${name}_sorted_markdupremove_sort_filter.bed

# stats the reads length distribution
awk '{print $3-$2}' ${name}_sorted_markdupremove_sort_filter.bed | sort -k1,1n |uniq -c | sed 's/\s\s*/ /g' | awk '{print $2"\t"$1}' > ${name}_readlengthcount_sorted_markdupremove_sort_filter.txt

# remove reads longer than 32nt
awk '{if($3-$2<=32) {print}}' ${name}_sorted_markdupremove_sort_filter.bed > ${name}_filter32_sorted_markdupremove_sort_filter.bed

# extract 26nt reads and get seq in bed
awk '{if($3-$2==26){print}}' ${name}_filter32_sorted_markdupremove_sort_filter.bed > ${name}_26nt_filter32_sorted_markdupremove_sort_filter.bed
bedtools getfasta -fi $ref \
	              -bed ${name}_26nt_filter32_sorted_markdupremove_sort_filter.bed \
	              -s \
	              -bedOut \
	              | awk '{ $7=toupper($7); OFS="\t"; print $0}' \
	              > ${name}_seq_26nt_filter32_sorted_markdupremove_sort_filter.bed

# remove reads containing N
awk -v OFS="\t" '$7!~/[^ATCG]/{print}' ${name}_seq_26nt_filter32_sorted_markdupremove_sort_filter.bed > ${name}_rmALL_seq_26nt_filter32_sorted_markdupremove_sort_filter.bed

# count total mapped reads
grep -c "^" ${name}_sorted_markdupremove_sort_filter.bed > ${name}_countreads_sorted_markdupremove_sort_filter.txt
grep -c "^" ${name}_filter32_sorted_markdupremove_sort_filter.bed > ${name}_readcount_filter32_sorted_markdupremove_sort_filter.txt

# bed to bdg
bedtools genomecov -i ${name}_filter32_sorted_markdupremove_sort_filter.bed \
                   -strand "+" \
                   -g $chrNamelength \
                   -bga \
                   | sort -k1,1 -k2,2n \
                   > ${name}_fwd_filter32_sorted_markdupremove_sort_filter.bdg

bedtools genomecov -i ${name}_filter32_sorted_markdupremove_sort_filter.bed \
                   -strand "-" \
                   -g $chrNamelength \
                   -bga \
                   | awk '{$4=-$4; OFS="\t"; print $0}' \
                   | sort -k1,1 -k2,2n \
                   > ${name}_rev_filter32_sorted_markdupremove_sort_filter.bdg

# normlize bdg by total mapping reads
NORMF=$(cat ${name}_filter32_sorted_markdupremove_sort_filter.bed | wc -l) # use 'cat' so the output will not contain file name
NORMF=$(python -c "print( $NORMF/(1000000*1000))") # RPKM=reads/(Million*kb)
awk -v NORMF=$NORMF '{$4=$4/NORMF; OFS="\t"; print $0}' ${name}_fwd_filter32_sorted_markdupremove_sort_filter.bdg > ${name}_RPKM_fwd_filter32_sorted_markdupremove_sort_filter.bdg
awk -v NORMF=$NORMF '{$4=$4/NORMF; OFS="\t"; print $0}' ${name}_rev_filter32_sorted_markdupremove_sort_filter.bdg > ${name}_RPKM_rev_filter32_sorted_markdupremove_sort_filter.bdg

# bdg into bw
bedGraphToBigWig ${name}_RPKM_fwd_filter32_sorted_markdupremove_sort_filter.bdg $chrNamelength ${name}_RPKM_fwd_filter32_sorted_markdupremove_sort_filter.bw
bedGraphToBigWig ${name}_RPKM_rev_filter32_sorted_markdupremove_sort_filter.bdg $chrNamelength ${name}_RPKM_rev_filter32_sorted_markdupremove_sort_filter.bw

# delete useless files
rm ${name}_filter.fastq
rm ${name}_filter.sai
rm ${name}_filter.sam
rm ${name}_filterpolyA.fastq
rm ${name}_filter.bam
rm ${name}_sort_filter.bam
rm ${name}_markdupremove_sort_filter.bam
rm ${name}_markdupremove_sort_filter.bed
rm ${name}_sorted_markdupremove_sort_filter.bed
#rm ${name}_filter32_sorted_markdupremove_sort_filter.bed
rm ${name}_26nt_filter32_sorted_markdupremove_sort_filter.bed
rm ${name}_seq_26nt_filter32_sorted_markdupremove_sort_filter.bed
rm *.bdg
rm *.bai

#damage.bed into bam, sort, build index
bedtools bedtobam -i ${name}_filter32_sorted_markdupremove_sort_filter.bed \
                  -g $chrNamelength \
                  > ${name}_filter32_sorted_markdupremove_sort_filter.bam
sambamba sort -t $cores -o ${name}_filter32_sorted_markdupremove_sort_filter.sorted.bam ${name}_filter32_sorted_markdupremove_sort_filter.bam
rm ${name}_filter32_sorted_markdupremove_sort_filter.bam
