#!/usr/bin/env bash

#SBATCH -J cisplatin_Damage-seq   
#SBATCH -p dna             
#SBATCH -N 1                 
#SBATCH --mem=80G                 
#SBATCH --cpus-per-task=20   
#SBATCH -t 03-00:00:00        
#SBATCH -o cisplatin_Damage-seq.Info.out      
#SBATCH --get-user-env

#--------使用说明--------
echo "该脚本使用命令需在测序数据文件夹Sample_xxxx所在目录输入,用法例下:"
echo "sbatch ~/Tanyuanqing/scripts/cisplatin_Damage-seq_scripts/slurm_auto_analysis_for_cisplatin_Damage-seq.sh human"
#-----------------------

#输入物种名称,确定分析脚本
if [ $1 = "human" ]
then
    script=~/Tanyuanqing/scripts/cisplatin_Damage-seq_scripts/cisplatin_Damage-seq_analysis_for_human_hg38.sh
elif [ $1 = "mouse" ]
then 
    script=~/Tanyuanqing/scripts/cisplatin_Damage-seq_scripts/cisplatin_Damage-seq_analysis_for_mouse_mm10.sh
else
    echo "物种名称参数缺失或输入错误，请重新运行脚本，并输入正确物种名称[human/mouse]";exit
fi

#提取测序数据文件夹名
ls -d Sample* > sn.list

#获取绝对路径
path=$PWD

#遍历每个数据文件夹里的read1和read2并进行脚本分析
for sl in $(cat sn.list)
do
     cd $path/$sl/ 
     read1=$(ls *R1.fastq.gz)
     read2=$(ls *R2.fastq.gz)
     bash $script $read1 $read2 20
     txt_2nt=$(ls *_count_seq_damageposition_R1_markdupremove_sort_fliter_discardAd1reads.txt)
     Rscript ~/Tanyuanqing/scripts/cisplatin_Damage-seq_scripts/cisplatindamage_2nt_precentage.plot.R $txt_2nt
done
