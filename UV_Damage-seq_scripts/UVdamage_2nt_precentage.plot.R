#!/opt/app/r/4.2.2/bin/Rscript
#Damage-seq位置2nt绘图

#通用设置
rm(list=ls())
options(stringsAsFactors = F)
#读取数据
args <- commandArgs(trailingOnly = TRUE)
cpd <- read.table(args[1])

#计算百分比
cpd["Ratio"] <- cpd$V2 / sum(cpd$V2) * 100
cpd <- cpd[,-2]
row <- nrow(cpd)

#绘图数据准备
library(reshape2)
cpd <- melt(cpd)
cpd$V1 <- factor(cpd$V1,levels = rev(cpd$V1))

#获取figure title
fig_title <- strsplit(args[1],split = "_")[[1]][1]

#ggplot绘图
library(ggplot2)
library(ggthemes)
colors <- c('#8dd3c7','#ffffb3','#bebada','#fb8072','#80b1d3',
            '#fdb462','#b3de69','#fccde5','#d9d9d9','#bc80bd',
            '#ccebc5','#ffed6f','#a6cee3','#1f78b4','#b2df8a',
            '#33a02c','#fb9a99','#e31a1c','#fdbf6f','#ff7f00',
            '#cab2d6','#6a3d9a','#ffff99','#b15928',"#E41A1C",
            "#1E90FF","#FF8C00","#4DAF4A","#984EA3","#40E0D0",
            "#FFC0CB","#00BFFF","#FFDEAD","#90EE90",
            "#EE82EE","#00FFFF","#F0A3FF", "#0075DC", 
            "#993F00","#4C005C","#2BCE48","#FFCC99",
            "#808080","#94FFB5","#8F7C00","#9DCC00",
            "#C20088","#003380","#FFA405","#FFA8BB",
            "#426600","#FF0010","#5EF1F2","#00998F",
            "#740AFF","#990000","#FFFF00")
p <- ggplot(data=cpd,aes(x=variable,y=value,fill=V1)) + geom_bar(stat = "identity",position = 'stack',width = 0.8) + 
     scale_y_continuous(labels = scales::percent_format(scale = 1)) +
     ggthemes::theme_base() +
     scale_fill_manual(values = sample(colors,row)) +
     ggtitle(fig_title) +
     theme(
           plot.title = element_text(size=6)
     ) +
     labs(
        y= "Dinucleotide frequency",
        x= NULL,
        fill= NULL
     )

#存储图片
ggsave(p,filename = paste(fig_title,".2nt.plot.pdf",sep=""),width = 3.8,device = cairo_pdf)
