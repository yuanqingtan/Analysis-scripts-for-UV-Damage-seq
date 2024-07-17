#!/usr/bin/env python3

import argparse
import os, sys
import pandas as pd
import matplotlib as mpl
mpl.use('Agg')
import pyBigWig as pbw
import multiprocessing as mp
from functools import reduce, partial
from collections import namedtuple


'''
All regions are 0-based.

usage:
plotProfile2.py gene.bed \
    -b sample1:sample1_mm10_F.bw,sample1_mm10_R.bw \
       sample2:sample2_mm10_F.bw,sample2_mm10_R.bw \
	--center region \
	--bins 100,300,100 \
	--lengths 5000,5000 \
	-o output

'''

# named tuples define
pairbw = namedtuple("pairbw", ["name", "f", "r"])
interval = namedtuple("interval", ["start","end","nbin"])
## regions will store interval
bedrecord = namedtuple("bedrecord", ["chrom","strand","name","regions"])

def get_args():
    arguments = argparse.ArgumentParser(description='''Plot profile from 
        stranded paired bigwig files.''')
    
    arguments.add_argument('bed', help='''At least 6-column tab-seperated
        bed file.''')

    arguments.add_argument('-o','--output', dest='out', 
        help='Output prefix.')
    
    arguments.add_argument('-b','--bigwig', nargs='+', 
        required=True, dest='bws', 
        help='''Bigwig pairs in 'Name:F.bw,R.bw' pattern.
            F.bw should have positive values and R.bw should have 
            negative values. eg. a:a_F.bw,a_R.bw b:b_F.bw,b_R.bw''')
    
    arguments.add_argument('--center', dest='center', default="start",
        choices=['start','end','region'],
        help='''Use which part in the bed file as the center.
                Default start.''')

    arguments.add_argument('--bins', dest='nbins', default="100,100",
        help='''Number of bins using for sum for each part. 
            If center is start and end, 2 numbers needed, eg. 100,200;
            If center is region, 3 numbers needed, eg 100,200,100.
            Default 100,100''')
    
    arguments.add_argument('--lengths', dest='lengths', 
        default="1000,1000",
        help='''Expanded length(bp) of each part to statistic.
            Two numbers needed, one for upstream, one for downstream,
            eg. 1000,2000. Default 1000,1000''')
    
    args = arguments.parse_args()

    # makdir 
    global mat_dir
    mat_dir = args.out+"_matrix"
    os.mkdir(mat_dir)
    # parse bin numbers
    bins = [int(i) for i in args.nbins.split(",")]
    args.up_bin = bins[0]
    args.down_bin = bins[-1]
    if len(bins) == 3:
        args.mid_bin = bins[1]
    
    # parse lengths
    lengths = [int(i) for i in args.lengths.split(",")]
    args.up_len = lengths[0]
    args.down_len = lengths[1]

    # parse bigwig pairs
    bw_names = []
    for i in args.bws:
        i = i.strip().split(":")
        ipair = i[1].split(',')
        bw_names.append( pairbw(i[0], ipair[0], ipair[1] ) )
    args.bws = bw_names

    # fetch chrom length
    def merge(x,y):
        d={}
        for k in ( x.keys() & y.keys() ):
            d[k] = min(x[k], y[k])
        return d
    
    chroms_all = []
    for i in args.bws:
        chroms_all.append( pbw.open(i.f).chroms() )
        chroms_all.append( pbw.open(i.r).chroms() )
    chromLength = reduce(merge, chroms_all)
    args.chroms = chromLength

    return args


def parse_bed( args ):

    '''
    Prepare bed file for summarizing.
    Final result is a list of namedtuple of bedrecord.
    argumments:
        args       args namespace
    '''

    bed_list = []

    def check( ainterval, chr_len):
        if (ainterval.start >0) and (ainterval.end <= chr_len) :
            return True
        else:
            return False  

    with open(args.bed) as bed:
        for line in bed:
            line = line.strip().split("\t")
            chrom = line[0]
            start = int(line[1])
            end = int(line[2])
            name = line[3]
            strand = line[5]

            if chrom not in args.chroms:
                continue

            if args.center == "start":
                if strand == '-':
                    up_interval = interval(end, end+args.up_len, args.up_bin)
                    down_interval = interval(end-args.down_len, end, args.down_bin)
                else:
                    up_interval = interval(start-args.up_len, start, args.up_bin)
                    down_interval = interval(start, start+args.down_len, args.down_bin)
                all_intervals = [up_interval, down_interval]
            elif args.center == "end":
                if strand == '-':
                    up_interval = interval(start, start+args.up_len, args.up_bin)
                    down_interval = interval(start-args.down_len, start, args.down_bin)
                else:
                    up_interval = interval(end-args.up_len, end, args.up_bin)
                    down_interval = interval(end, end+args.down_len, args.down_bin)
                all_intervals = [up_interval, down_interval]
            elif args.center == "region":
                mid_interval = interval(start, end, args.mid_bin)
                if strand == '-':
                    up_interval = interval(end, end+args.up_len, args.up_bin)
                    down_interval = interval(start-args.down_len, start, args.down_bin)
                else:
                    up_interval = interval(start-args.up_len, start, args.up_bin)
                    down_interval = interval(end, end+args.down_len, args.down_bin)
                all_intervals = [up_interval, mid_interval, down_interval]
            
            abedline = bedrecord(chrom, strand, name, all_intervals)
            chr_len = args.chroms[chrom]
            check_res = [check(i, chr_len) for i in all_intervals]
            if all(check_res):
                bed_list.append(abedline)

    return(bed_list)


def count( bw, bed ):
    
    bw = pairbw( bw.name, pbw.open(bw.f), pbw.open(bw.r) )

    f_matrix = []
    r_matrix = []
    for ibed in bed:
        i_f = [ibed.name]
        i_r = [ibed.name]
        if ibed.strand == "-":
            for r in ibed.regions:
                i_f += bw.r.stats(ibed.chrom, r.start, r.end, nBins=r.nbin, exact=True)[::-1]
                i_r += bw.f.stats(ibed.chrom, r.start, r.end, nBins=r.nbin, exact=True)[::-1]
        else:
            for r in ibed.regions:
                i_f += bw.f.stats(ibed.chrom, r.start, r.end, nBins=r.nbin, exact=True) 
                i_r += bw.r.stats(ibed.chrom, r.start, r.end, nBins=r.nbin, exact=True)
        
        f_matrix.append(i_f)
        r_matrix.append(i_r)

    # write 
    f_matrix = pd.DataFrame(f_matrix).fillna(0).set_index(0).abs()
    f_matrix.to_csv(f'{mat_dir}/{bw.name}_F.matrix.csv', header=False)
    r_matrix = -pd.DataFrame(r_matrix).fillna(0).set_index(0).abs()
    r_matrix.to_csv(f'{mat_dir}/{bw.name}_R.matrix.csv', header=False)

    # mean
    ave_dict = {bw.name+"_F": f_matrix.mean(0),
                bw.name+"_R": r_matrix.mean(0)}
    
    return ave_dict

def main():
    args = get_args()
    bed_list = parse_bed( args )

    threads = mp.Pool()
    count_bed = partial(count, bed=bed_list)

    ave_dict = {}
    for a_dict in threads.imap( count_bed, args.bws):
        ave_dict.update(a_dict)

    ave_df = pd.DataFrame(ave_dict)
    ave_df.to_csv(args.out+".csv", index=False)
    ave_ax = ave_df.plot.line()
    ave_ax.figure.savefig(args.out+".pdf")

if __name__ == "__main__":
    main()
