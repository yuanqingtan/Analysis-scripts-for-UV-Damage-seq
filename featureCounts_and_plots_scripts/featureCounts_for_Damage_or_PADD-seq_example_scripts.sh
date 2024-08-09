#!/bin/bash

gtf_hg38=~/Tanyuanqing/GTF_GFF3files/hg38GTF_GFF3/gencode.v40.annotation.gtf

featureCounts -a $gtf_hg38 \
              -o gene.readscount.TS.tab ELOF1-KO.0h.bam ELOF1-KO.8h.bam UVSSA-KO.0h.bam UVSSA-KO.8h.bam \
                 UVSSA-KO-WT.0h.bam UVSSA-KO-WT.8h.bam UVSSA-KO-K414R.0h.bam UVSSA-KO-K414R.8h.bam \
              -t gene \
              -T 8 \
              -s 2 \
              --minOverlap 2 \
              -g gene_name

featureCounts -a $gtf_hg38 \
              -o gene.readscount.NTS.tab ELOF1-KO.0h.bam ELOF1-KO.8h.bam UVSSA-KO.0h.bam UVSSA-KO.8h.bam \
                 UVSSA-KO-WT.0h.bam UVSSA-KO-WT.8h.bam UVSSA-KO-K414R.0h.bam UVSSA-KO-K414R.8h.bam \
              -t gene \
              -T 8 \
              -s 1 \
              --minOverlap 2 \
              -g gene_name
              
