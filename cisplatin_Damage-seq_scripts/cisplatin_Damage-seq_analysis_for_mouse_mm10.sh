#!/usr/bin/env bash    

#---------------------------------------------------------------------------
# Version:   24.07.01
# Author:    Tan yuanqing
# Function:  cisplatin Damage-seq analysis for mouse mm10
# Usage:     bash cisplatin_Damage-seq_analysis_for_mouse_mm10.sh Read1 Read2 cores
#---------------------------------------------------------------------------   

# 所需参数和变量设置
name=$(basename $1 "_R1.fastq.gz")
name1=$(basename $1 ".fastq.gz")
name2=$(basename $2 ".fastq.gz")
chrNamelength=~/Tanyuanqing/Genome/mm10/bwa_index/mm10chrNamelength.txt
ref=~/Tanyuanqing/Genome/mm10/bwa_index/mm10.fa
cores=$3


# 去除文库中的非特异性reads
cutadapt -g GACTGGTTCCAATTGAAAGTGCTCTTCCGATCT \
         --discard-trimmed \
         --pair-filter=any \
         -j $cores \
         -o ${name1}_discardAd1reads.fastq.gz \
         -p ${name2}_discardAd1reads.fastq.gz \
         $1 $2

# 去除adapt
fastp -i ${name1}_discardAd1reads.fastq.gz \
      -I ${name2}_discardAd1reads.fastq.gz \
      -o ${name1}_fliter_discardAd1reads.fastq.gz \
      -O ${name2}_fliter_discardAd1reads.fastq.gz \
      -w $cores \
      -q 30 

# mapping by BWA
bwa mem -t $cores \
           $ref \
           ${name1}_fliter_discardAd1reads.fastq.gz \
           ${name2}_fliter_discardAd1reads.fastq.gz \
           > ${name}_fliter_discardAd1reads.sam

# sam into bam, sort and remove duplicated
samtools view -b -S -@ $cores ${name}_fliter_discardAd1reads.sam > ${name}_fliter_discardAd1reads.bam
sambamba sort -t $cores -o ${name}_sort_fliter_discardAd1reads.bam ${name}_fliter_discardAd1reads.bam
sambamba markdup -r -t $cores ${name}_sort_fliter_discardAd1reads.bam ${name}_markdupremove_sort_fliter_discardAd1reads.bam

# 提取read1，并过滤低mapping质量的reads
samtools view -b -f 0x40 -q 25 -o ${name}_R1_markdupremove_sort_fliter_discardAd1reads.bam ${name}_markdupremove_sort_fliter_discardAd1reads.bam

# bam into bed
bedtools bamtobed -i ${name}_R1_markdupremove_sort_fliter_discardAd1reads.bam > ${name}_R1_markdupremove_sort_fliter_discardAd1reads.bed

# 提取损伤坐标
bedtools flank -i ${name}_R1_markdupremove_sort_fliter_discardAd1reads.bed \
               -g $chrNamelength \
               -l 2 -r 0 -s > ${name}_damageposition_R1_markdupremove_sort_fliter_discardAd1reads.bed

# 提取损伤位置碱基序列
bedtools getfasta -fi $ref \
                  -bed ${name}_damageposition_R1_markdupremove_sort_fliter_discardAd1reads.bed \
                  -s \
                  -bedOut \
                  | awk '{$7=toupper($7); OFS="\t"; print $0}' \
                  > ${name}_seq_damageposition_R1_markdupremove_sort_fliter_discardAd1reads.bed

# 计算损伤位点二嘧啶的数量和比例
cut -f7 ${name}_seq_damageposition_R1_markdupremove_sort_fliter_discardAd1reads.bed|sort|uniq -c|sort -k1r,1n \
    | awk '{OFS="\t";print $2,$1}' > ${name}_count_seq_damageposition_R1_markdupremove_sort_fliter_discardAd1reads.txt

# 过滤得到特异性的损伤reads
awk '$7 ~ /GG/' ${name}_seq_damageposition_R1_markdupremove_sort_fliter_discardAd1reads.bed | \
    sort -k1,1 -k2,2n > ${name}_cisplatin_seq_damageposition_R1_markdupremove_sort_fliter_discardAd1reads.bed

# 将正链reads和负链reads分开计算各自的coverage，生成相应的bdg文件
bedtools genomecov -i ${name}_cisplatin_seq_damageposition_R1_markdupremove_sort_fliter_discardAd1reads.bed \
                   -strand "+" \
                   -g $chrNamelength \
                   -bga \
                   | sort -k1,1 -k2,2n \
                   > ${name}_fwd_cisplatin_seq_damageposition_R1_markdupremove_sort_fliter_discardAd1reads.bdg 

bedtools genomecov -i ${name}_cisplatin_seq_damageposition_R1_markdupremove_sort_fliter_discardAd1reads.bed \
                   -strand "-" \
                   -g $chrNamelength \
                   -bga \
                   | awk '{$4=-$4; OFS="\t"; print $0}' \
                   | sort -k1,1 -k2,2n \
                   > ${name}_rev_cisplatin_seq_damageposition_R1_markdupremove_sort_fliter_discardAd1reads.bdg

# normlize coverage value with RPKM
NORMF=$(cat ${name}_cisplatin_seq_damageposition_R1_markdupremove_sort_fliter_discardAd1reads.bed |wc -l) # use 'cat' so the output will not contain file name
NORMF=$(python -c "print( $NORMF/(1000000*1000))") # RPKM=reads/(Million*kb)
awk -v NORMF=$NORMF '{$4=$4/NORMF; OFS="\t"; print $0}' ${name}_fwd_cisplatin_seq_damageposition_R1_markdupremove_sort_fliter_discardAd1reads.bdg > ${name}_normalize_fwd_cisplatin_seq_damageposition_R1_markdupremove_sort_fliter_discardAd1reads.bdg
awk -v NORMF=$NORMF '{$4=$4/NORMF; OFS="\t"; print $0}' ${name}_rev_cisplatin_seq_damageposition_R1_markdupremove_sort_fliter_discardAd1reads.bdg > ${name}_normalize_rev_cisplatin_seq_damageposition_R1_markdupremove_sort_fliter_discardAd1reads.bdg

# bdg into bigwig
bedGraphToBigWig ${name}_normalize_fwd_cisplatin_seq_damageposition_R1_markdupremove_sort_fliter_discardAd1reads.bdg $chrNamelength ${name}_RPKM_normalize_fwd_cisplatin_seq_damageposition_R1_markdupremove_sort_fliter_discardAd1reads.bw
bedGraphToBigWig ${name}_normalize_rev_cisplatin_seq_damageposition_R1_markdupremove_sort_fliter_discardAd1reads.bdg $chrNamelength ${name}_RPKM_normalize_rev_cisplatin_seq_damageposition_R1_markdupremove_sort_fliter_discardAd1reads.bw

# count total mapped reads
wc -l ${name}_R1_markdupremove_sort_fliter_discardAd1reads.bed > ${name}_readscountTotal_R1_markdupremove_sort_fliter_discardAd1reads.txt

# delete useless files
rm ${name1}_discardAd1reads.fastq.gz
rm ${name2}_discardAd1reads.fastq.gz
rm ${name1}_fliter_discardAd1reads.fastq.gz
rm ${name2}_fliter_discardAd1reads.fastq.gz
rm ${name}_fliter_discardAd1reads.sam
rm *fastp*
rm *.bdg
#rm *.bed
rm ${name}_seq_damageposition_R1_markdupremove_sort_fliter_discardAd1reads.bed
rm ${name}_damageposition_R1_markdupremove_sort_fliter_discardAd1reads.bed
rm ${name}_R1_markdupremove_sort_fliter_discardAd1reads.bed
rm *.bam
rm *.bai

# damage.bed into bam, sort and build index
bedtools bedtobam -i ${name}_cisplatin_seq_damageposition_R1_markdupremove_sort_fliter_discardAd1reads.bed \
                  -g $chrNamelength \
                  > ${name}_cisplatin_seq_damageposition_R1_markdupremove_sort_fliter_discardAd1reads.bam
sambamba sort -t $cores -o ${name}_cisplatin_seq_damageposition_R1_markdupremove_sort_fliter_discardAd1reads.sorted.bam \
                           ${name}_cisplatin_seq_damageposition_R1_markdupremove_sort_fliter_discardAd1reads.bam
rm ${name}_cisplatin_seq_damageposition_R1_markdupremove_sort_fliter_discardAd1reads.bam
