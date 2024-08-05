#!/usr/bin/env python
# coding: utf-8

# # XRseq reads 连续2nt 绘图

# ###导入模块

# In[6]:


import matplotlib.pyplot as plt
import matplotlib as mpl
import pandas as pd
import numpy as np
import os
import sys
#import seaborn as sns


# #设置工作文件夹

# In[7]:

# #获取当前文件夹csv文件名list

# In[8]:


#创建绘图颜色list
color_list = ['#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#ffff99','#b15928','#b3e2cd','#fdcdac','#cbd5e8','#f4cae4']


# #定义数据读取及绘图函数

# In[9]:


def rdata_plot(sn):
    
    df = pd.read_csv(sn,index_col = 0)*100#读取数据
    x = range(len(df.index))#获取x轴数据
    c_list = df.columns.to_list()#获取column name list

    #设置输出文件名

    main = sn.split('.')[:-1]
    main = ''.join(main)
    #开始绘图设置各种参数
    fig, ax = plt.subplots(dpi=100,figsize=(10,6))
    #全局设置字体大小粗细
    plt.rc('font',family='Arial',weight = 'bold',size=12)
    ##给X轴生成文本标签
    xtick = []
    for i in range(0,len(df.index),2):
        xtick.append(str(i+1) + '-' + str(i+2))
    bottom = 0
    ##创建堆叠图
    for i in range(len(c_list)):
        ax.bar(x,df[c_list[i]],label = c_list[i],bottom = bottom,width=1,edgecolor="white", linewidth=0.5,color=color_list[i])
        bottom = bottom + df[c_list[i]]
    ax.legend()
    ax.spines['left'].set_color(None)
    ax.spines[['top','right']].set_visible(False)
    ax.spines['bottom'].set_position(('data',-1))
    ax.set_xlim(-0.55,len(df.index)-0.5)
    ax.set_title(label=main,fontsize=16,color='black',weight='bold')
    ax.set_ylabel('Dinucleotide frequency (%)',fontsize=16,color='black',weight='bold')
    ax.set_xlabel('Position',fontsize=16,color='black',weight='bold')
    ax.set_xticks(range(0,len(df.index),2),xtick,rotation = 30)
    ax.legend(bbox_to_anchor=(1.01, 0), loc='lower left', borderaxespad=0,prop ={'size':13},
              ncol = 1,frameon=False,labelcolor=color_list)
    
    #输出图片
    plt.savefig(main +'.pdf')


# In[10]:
if __name__ == "__main__":
    rdata_plot(sys.argv[1])

