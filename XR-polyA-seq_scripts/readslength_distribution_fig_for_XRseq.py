#!/usr/bin/env python
# coding: utf-8

# # reads length distrbution figure

# In[11]:


import matplotlib.pyplot as plt
import matplotlib as mpl
import pandas as pd
import numpy as np
import os
import sys


#定义数据读取绘图函数
def rdata_plot(sn):
##读取文件

    df = pd.read_table(sn,sep='\t',header=None)
    ##设置输出文件名
    main = sn.split('.')[:-1]
    main = ''.join(main)
    ##获取x轴，y轴数据

    x = df.iloc[:30,0]
    #y_sum = df.iloc[:,1]
    y = df.iloc[:30,1]
    y = (y/sum(y))*100

    ##绘图

    fig, ax = plt.subplots(figsize=(8,6),dpi=150)##创建一个含有单一坐标轴的图
    #全局设置字体大小粗细
    plt.rc('font',family='Arial',size=12)
    ax.tick_params(labelsize=10)
    ax.bar(x,y,width=1, edgecolor="black", linewidth=0.5,color = '#8dd3c7',alpha=1)
    ax.set_title(label=main,fontsize=18,color='black')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_position(('data',-0.25))
    ax.set_ylim(0,max(y))
    ax.set_ylabel('Read count (%)',fontsize=14)
    ax.set_yticks([0,5,10,15])
    ax.spines['left'].set_color(None)
    ax.set_xlabel('Read length',fontsize=14)
    ax.set_xlim(min(x)-1,max(x)+1)
    ax.set_xticks([15,20,25,30,35,40])
    ##输出图片文件
    plt.savefig(main + '.pdf')


# In[15]:
if __name__ == "__main__":
    rdata_plot(sys.argv[1])

