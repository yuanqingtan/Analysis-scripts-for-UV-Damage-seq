#!/usr/bin/env bash 

#SBATCH -J plotprofile2   
#SBATCH -p interactive             
#SBATCH -N 1                 
#SBATCH --mem=5G                 
#SBATCH --cpus-per-task=1   
#SBATCH -t 0-06:00:00        
#SBATCH -o Info.plotprofile2.out      
#SBATCH --get-user-env


# 获取样本数据文件夹
ls -d Sample* > sn.list

#根据输入参数选择bed文件

case $1 in 
all) bed=~/Tanyuanqing/GTF_GFF3files/hg38GTF_GFF3/changeorder_gencode.v40.annotation_GeneProteinCodingLevel12_nonoverlapped_sortedbygenelength.bed;;
high) bed=~/Tanyuanqing/GTF_GFF3files/hg38GTF_GFF3/Helacellrnaseq/changeorder_Hela_gencode.v40.annotation_GeneProteinCodingLevel12_nonoverlapped_5+TPM_3901_sortedbygenelength.bed;;
low) bed=~/Tanyuanqing/GTF_GFF3files/hg38GTF_GFF3/Helacellrnaseq/changeorder_Hela_gencode.v40.annotation_GeneProteinCodingLevel12_nonoverlapped_0-5TPM_3571_sortedbygenelength.bed;;
no) bed=~/Tanyuanqing/GTF_GFF3files/hg38GTF_GFF3/Helacellrnaseq/changeorder_Hela_gencode.v40.annotation_GeneProteinCodingLevel12_nonoverlapped_0TPM_2287_sortedbygenelength.bed;;
allhe) bed=~/Tanyuanqing/GTF_GFF3files/hg38GTF_GFF3/Helacellrnaseq/PADD-seq.proteincoding_highexpressed.nonoverlapped.bed/changeorder.gencode.v40.annotation_GeneProteinCodingLevel12.HelaTPM5+.nonoverlapped.sortedbygenelength.5423.bed;;
*) echo "请输入正确的参数[all,high,low,no]";exit;;
esac

#设置plotprofile2.py脚本所需变量
bed_type=${1}bed      #bed文件类型描述,可选[allbed,highbed,lowbed,nonbed]
script=~/Tanyuanqing/scripts/plotProfile2.py #脚本路径
center="region"       #可选["start","region","end"],依次代表TSS,genebody,TTS
bins="5,50,5"    #设置区段均分数
lengths="2000,2000"   #设置center两侧的区段范围

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
    cd  $path/$sl/
    bw_F=$(ls *fwd*.bw)
    bw_R=$(ls *rev*.bw)
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
    
