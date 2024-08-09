#!/usr/bin/env bash 

#SBATCH -J plotprofile2   
#SBATCH -p interactive             
#SBATCH -N 1                 
#SBATCH --mem=5G                 
#SBATCH --cpus-per-task=1   
#SBATCH -t 0-06:00:00        
#SBATCH -o Info.plotprofile2.out      
#SBATCH --get-user-env


#sample fold list
ls -d Sample* > sn.list

#根据输入参数选择bed文件

case $1 in 
all) bed=~/Tanyuanqing/GTF_GFF3files/hg38GTF_GFF3/changeorder_gencode.v40.annotation_GeneProteinCodingLevel12_nonoverlapped_sortedbygenelength.bed;;
high) bed=~/Tanyuanqing/GTF_GFF3files/hg38GTF_GFF3/XP-C_RNAseq/changeorder_gencode.v40.annotation_GeneProteinCodingLevel12_nonoverlapped_1+TPM_3039_sortedbygenelength.bed;;
low) bed=~/Tanyuanqing/GTF_GFF3files/hg38GTF_GFF3/XP-C_RNAseq/changeorder_gencode.v40.annotation_GeneProteinCodingLevel12_nonoverlapped_0-1TPM_3245_sortedbygenelength.bed;;
no) bed=~/Tanyuanqing/GTF_GFF3files/hg38GTF_GFF3/XP-C_RNAseq/changeorder_gencode.v40.annotation_GeneProteinCodingLevel12_nonoverlapped_0TPM_3489_sortedbygenelength.bed;;
alle) bed=~/Tanyuanqing/GTF_GFF3files/hg38GTF_GFF3/XP-C_RNAseq/XP-C_PADD-seq_genebed/changeorder.gencode.v40.annotation_GeneProteinCodingLevel12.0+TPM13841.nonoverlapped7318.sortedbygenelength.bed;;
allhe) bed=~/Tanyuanqing/GTF_GFF3files/hg38GTF_GFF3/XP-C_RNAseq/XP-C_PADD-seq_genebed/changeorder.gencode.v40.annotation_GeneProteinCodingLevel12.1+TPM.nonoverlapped4869.sortedbygenelength.bed;;
he40k+) bed=~/Tanyuanqing/GTF_GFF3files/hg38GTF_GFF3/XP-C_RNAseq/XP-C_PADD-seq_genebed/changeorder.gencode.v40.annotation_GeneProteinCodingLevel12.1+TPM6946.40k+3328.nonoverlapped2881.sortedbygenelength.bed;;
e40k+) bed=~/Tanyuanqing/GTF_GFF3files/hg38GTF_GFF3/XP-C_RNAseq/XP-C_PADD-seq_genebed/changeorder.gencode.v40.annotation_GeneProteinCodingLevel12.0+TPM13841.40k+6309.nonoverlapped4816.sortedbygenelength.bed;;
TPM1+) bed=~/Tanyuanqing/GTF_GFF3files/hg38GTF_GFF3/XP-C_RNAseq/XP-C_PADD-seq_genebed/gencode.v40.annotation_GeneProteinCodingLevel12.1+TPM6946.bed;;
TPM0+) bed=~/Tanyuanqing/GTF_GFF3files/hg38GTF_GFF3/XP-C_RNAseq/XP-C_PADD-seq_genebed/gencode.v40.annotation_GeneProteinCodingLevel12.0+TPM13841.bed;;
allhee) bed=~/Tanyuanqing/GTF_GFF3files/hg38GTF_GFF3/XP-C_RNAseq/XP-C_PADD-seq_genebed/changeorder.gencode.v40.annotation_GeneProteinCodingLevel12.1+TPM.nonoverlapped4869.bed;;
*) echo "请输入正确的参数[all,high,low,no] for Damage-seq, [alle,allhe,he40k+,e40k+] for PADD-seq,选择目标gene bed 文件进行计算";exit;;
esac

#设置plotprofile2.py脚本所需变量
bed_type=${1}bed      #bed文件类型描述,可选[allbed,highbed,lowbed,nonbed]
script=~/Tanyuanqing/scripts/plotProfile2.py #脚本路径
center="region"       #可选["start","region","end"],依次代表TSS,genebody,TTS
bins="10,30,10"    #设置区段均分数
lengths="5000,5000"   #设置center两侧的区段范围

#打印所选参数信息
echo "打印参数信息如下："
echo "script=$script"
echo "bed=$1:$bed"
echo "center=$center"
echo "bins=$bins"
echo "lengths=$lengths"

#设置数据助记变量
if [ $center = "region" ]
then
    refpoint="TSS-TTS"
elif [ $center = "start" ]
then
    refpoint="TSS"
elif [ $center = "end" ]
then
    refpoint="TTS"
else
    echo "参数center选择错误，不在可选范围["start","region","end"]内";exit
fi

Bins=`echo $bins | sed 's/,/_/g'`
Bins=Bins${Bins}
Lengths=`echo $lengths | sed 's/,/_/g'`
Lengths=lengths${Lengths}

path=$PWD #获取当前绝对路径

#依次进入各个数据文件夹,寻找到bw文件,运行脚本

for sl in $(cat sn.list)
do
    cd $path/$sl/
    bw_F=$(ls *fwd*reads.bw)
    bw_R=$(ls *rev*reads.bw)
    $script $bed \
            -b sample1:$bw_F,$bw_R \
	    --center $center \
	    --bins $bins \
	    --lengths $lengths \
	    -o ${sl}_${bed_type}_${refpoint}_${Bins}_${Lengths}_plotprofile2
    mv ${sl}_${bed_type}_${refpoint}_${Bins}_${Lengths}_plotprofile2.* ${sl}_${bed_type}_${refpoint}_${Bins}_${Lengths}_plotprofile2_matrix/
    cd ${sl}_${bed_type}_${refpoint}_${Bins}_${Lengths}_plotprofile2_matrix/
    mv sample1_F.matrix.csv ${sl}_${bed_type}_${refpoint}_${Bins}_${Lengths}_sample1_F.matrix.csv
    mv sample1_R.matrix.csv ${sl}_${bed_type}_${refpoint}_${Bins}_${Lengths}_sample1_R.matrix.csv
done
    
