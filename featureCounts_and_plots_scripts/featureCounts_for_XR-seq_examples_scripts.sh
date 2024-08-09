#!/bin/bash

gtf_hg38=~/Tanyuanqing/GTF_GFF3files/hg38GTF_GFF3/gencode.v40.annotation.gtf

featureCounts -a $gtf_hg38 \
              -o gene.readscount.TS.tab WT.rep1.bam WT.rep2.bam CSBKO.bam ELOF1KO.bam STK19KO32.rep1.bam STK19KO48.rep1.bam  \
                 STK19KO32.rep2.bam STK19KO48.rep2.bam STK19KO32-C1.bam STK19KO32-C2.bam STK19KO48-C1.bam STK19KO48-C2.bam \
                 STK19KO32-C1-K317P.bam STK19KO32-C2-K317P.bam STK19KO48-C1-K317P.bam STK19KO48-C2-K317P.bam \
              -t gene \
              -T 8 \
              -s 2 \
              --minOverlap 14 \
              -g gene_name

featureCounts -a $gtf_hg38 \
              -o gene.readscount.NTS.tab WT.rep1.bam WT.rep2.bam CSBKO.bam ELOF1KO.bam STK19KO32.rep1.bam STK19KO48.rep1.bam  \
                 STK19KO32.rep2.bam STK19KO48.rep2.bam STK19KO32-C1.bam STK19KO32-C2.bam STK19KO48-C1.bam STK19KO48-C2.bam \
                 STK19KO32-C1-K317P.bam STK19KO32-C2-K317P.bam STK19KO48-C1-K317P.bam STK19KO48-C2-K317P.bam \
              -t gene \
              -T 8 \
              -s 1 \
              --minOverlap 14 \
              -g gene_name
              
