#!/bin/bash

bam=$1
bz=$2
sml=$3

bamCoverage -b $bam \
            -o ${bam}_binsize${bz}_smoothLength${sml}.fwd.RPKM.bw \
            --samFlagExclude 16 \
            --binSize $bz \
            --smoothLength $sml \
            -p 20 \
            --normalizeUsing RPKM
            

bamCoverage -b $bam \
            -o ${bam}_binsize${bz}_smoothLength${sml}.rev.RPKM.bw \
            --samFlagInclude 16 \
            --scaleFactor -1 \
            --binSize $bz \
            --smoothLength $sml \
            -p 20 \
            --normalizeUsing RPKM
