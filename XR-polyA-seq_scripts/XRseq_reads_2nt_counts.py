#!/usr/bin/env python
# coding: utf-8

# # XRseq reads 连续2nt统计

# #此处使用的bed文件要求第7列为reads序列，若不是，则需修改以下代码中df[7]中的数字
# 注意：此处reads是不含有碱基为N的reads，一般在生成含reads序列的bed文件之后，要过滤掉含N的reads

# ###导入模块

import pandas as pd
import numpy as np
import os
import sys

#定义数据读取统计函数
def rdata_count(sn):

###读取文件

    df = pd.read_table(sn,sep = '\t',header = None)
###设置输出文件名
    main = sn.split('.')[:-1]
    main = ''.join(main)
###取第7列reads序列

    df_seq = df[6]

###获取行数=reads数

    row_number = len(df_seq.index)

###将每行reads序列列表化，并生成一个DataFrame

    mdict = []
    for i in range(row_number):
        mdict.append(list(df_seq[i]))

    DF = pd.DataFrame(mdict)
###获取列数
    col_number = len(DF.columns)

###统计相邻位置不同类别2碱基个数

    c = []
    for i in range(col_number - 1):
        DF_2nt = DF[i] + DF[i+1]
        c.append(DF_2nt.value_counts([0]))

    c = pd.DataFrame(c,columns = ["TT","TC","CT","CC","AG","CA","AA","AC","AT","GC","GA","GT","TG","GG","TA","CG"])
    c.to_csv(f'{main}.csv')

if __name__ == "__main__":
    rdata_count(sys.argv[1])

