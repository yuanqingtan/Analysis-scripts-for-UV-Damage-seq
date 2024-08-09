#feature_counts.plot

#
rm(list=ls())
setwd("~/Tanyuanqing/project/Damage-seq/XP-C/XP-C_UV8h_0h/ELOF1_UVSSA_UVSSA-WT_UVSSA-K414R/s1608g01179_SHDWYY-20240422-L-02-2024-04-300939/featureCounts/")

#读取gene bed文件
genes <- read.table("gencode.v40.annotation_GeneProteinCodingLevel12.1+TPM6946.bed")

#读取数据
TS <- read.table("gene.readscount.TS.tab",header = T)
TS <- TS[TS$Geneid %in% genes$V4,]
TS <- TS[,-c(1:6)]

NTS <- read.table("gene.readscount.NTS.tab",header = T)
NTS <- NTS[NTS$Geneid %in% genes$V4,]
NTS <- NTS[,-c(1:6)]

div <- NTS/TS
div <- data.frame(div$ELOF1.KO.0h.bam/div$ELOF1.KO.8h.bam, div$UVSSA.KO.0h.bam/div$UVSSA.KO.8h.bam,
                  div$UVSSA.KO.WT.0h.bam/div$UVSSA.KO.WT.8h.bam, div$UVSSA.KO.K414R.0h.bam/div$UVSSA.KO.K414R.8h.bam)

div <- log2(div)
div <- na.omit(div)
div <- div[abs(rowSums(div))!=Inf,]
div <- na.omit(div)
names(div) <- c("ELOF1-KO","UVSSA-KO","UVSSA-KO-WT","UVSSA-KO-K414R")

# 假设df是你的dataframe，并且你想要检查所有列
# 找出所有行中至少有一个数值不在-1到1之间的行
rows_to_keep <- apply(div, 1, function(x) all(x >=-1 & x <= 1))

# 使用逻辑索引来选择这些行
div <- div[rows_to_keep, ]

# 现在div包含了所有没有数值在-1到1之间的行

#转换数据
library(reshape2)
library(ggpubr)
library(ggthemes)
library(patchwork)
# #low <- "#e5f5f9"
# low <- "#7FC97F"
# #high <- "#2ca25f"
# high <- "#BEAED4"
colors <- c("#fee8c8","#fdbb84","#efedf5","#bcbddc") 
div.re <- melt(div,value.name = 'log2(NTS/TS)(0h/8h)',variable.name = 'Sample')

#绘图
fig_title <- "high Genes featureCounts of Damage-seq in XP-C " #图标题
num <- nrow(div) 
annotation <- paste('(n = ',num,')')            #基因个数文本注释

p1 <- ggplot(div.re,aes(x=Sample,y=`log2(NTS/TS)(0h/8h)`,fill=Sample)) + geom_violin(trim = T,linewidth=0.2) +
  theme_classic() +
  geom_boxplot(width=0.2,cex=0.2,outlier.shape = NA) +
  geom_hline(aes(yintercept = 0),linetype = "dashed",linewidth=0.25,col="red",alpha=1) +
  ggtitle(fig_title) +
  stat_summary(fun = mean,geom = "point", shape=23,size = 2, color = "black",fill="white") +
  scale_fill_manual(values = colors) +
  #ggthemes::theme_base() +
  theme(axis.text.x = element_text(angle = 60,hjust = 1),
        plot.title = element_text(size = 10),
        #axis.ticks.x = element_blank(),
        panel.grid = element_blank(),
        
        ) +
  labs(
    y="log2(TS(8h/0h)/NTS(8h/0h))"
  ) +
  annotate("text", x = -Inf, y = Inf, size = 3.5, label = annotation,hjust = -.1,vjust = 1.5,color = "red") #+
  #scale_y_continuous(limits = c(-1,1))

# p1 <- p1 + stat_summary(fun = "mean", geom = "text", 
#                         aes(label = round(..y.., 4)), 
#                         col = "red", size = 2, fontface = "bold", vjust = 50)

my_c <- list(c("ELOF1-KO","UVSSA-KO"),c("UVSSA-KO","UVSSA-KO-WT"),c("UVSSA-KO-WT","UVSSA-KO-K414R"),
             c("UVSSA-KO","UVSSA-KO-K414R"),c("ELOF1-KO","UVSSA-KO-WT"),c("ELOF1-KO","UVSSA-KO-K414R"))
p1 <- p1+stat_summary(fun = "mean", geom = "text", 
                aes(label = round(..y.., 4)), 
                col = "red", size = 4, fontface = "bold", vjust = 5) +
   stat_compare_means(comparisons = my_c,method = "t.test",paired=T,label = "p.adj")
  
p1

# #test
pV1=t.test(div$`ELOF1-KO`,div$`UVSSA-KO`,paired=T)
# #format
pV1=format(pV1,digits = 3, scientific=T)
pV1="2.05e-04"

pV2=t.test(div$`UVSSA-KO`,div$`UVSSA-KO-WT`,paired = T)
pv2=format(pV2,digits = 3, scientific=T)
pv2="0"

pv3=t.test(div$`UVSSA-KO-WT`,div$`UVSSA-KO-K414R`,paired=T)
pv3=format(pv3,digits = 3, scientific=T)
pv3="3.52e-196"

pv4=t.test(div$`ELOF1-KO`,div$`UVSSA-KO-WT`,paired = T)
pv4=format(pv4,digits = 3, scientific=T)
pv4="0"

pv5=t.test(div$`ELOF1-KO`,div$`UVSSA-KO-K414R`,paired=T)
pv5=format(pv5,digits = 3, scientific=T)
pv5="1.03e-40"

pv6=t.test(div$`UVSSA-KO`,div$`UVSSA-KO-K414R`,paired = T)
pv6=format(pv6,digits = 3, scientific=T)
pv6="1.34-e25"

# W_C=5.107956e-243
# C_S=2.6e-02
# S_SC=1.2e-180
# W_SC=1e-34
# W_S=
# C_SC=

