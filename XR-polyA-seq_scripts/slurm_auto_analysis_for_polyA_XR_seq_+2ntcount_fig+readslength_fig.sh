#!/usr/bin/env bash

#SBATCH -J XR-seq-polyA   
#SBATCH -p rna             
#SBATCH -N 1                 
#SBATCH --mem=50G                 
#SBATCH --cpus-per-task=10  
#SBATCH -t 03-00:00:00        
#SBATCH -o XR-seq-polyA.Info.out      
#SBATCH --get-user-env

#--------使用说明--------
echo "该脚本使用命令需在测序数据文件夹Sample_xxxx所在目录输入,用法例下:"
echo "sbatch ~/Tanyuanqing/scripts/XR-polyA-seq_scripts/slurm_auto_analysis_for_polyA_XR_seq_+2ntcount_fig+readslength_fig.sh human"
#-----------------------

#输入物种名称,确定分析脚本
if [ $1 = "human" ]
then
    script=~/Tanyuanqing/scripts/XR-polyA-seq_scripts/XR-polyA-seq_analysis_for_human_hg38.sh
elif [ $1 = "mouse" ]
then 
    script=~/Tanyuanqing/scripts/XR-polyA-seq_scripts/XR-polyA-seq_analysis_for_mouse_mm10.sh
elif [ $1 = "dm6" ]
then
    script=~/Tanyuanqing/scripts/XR-polyA-seq_scripts/XR-polyA-seq_analysis_for_dm6.sh
else
    echo "物种名称输入错误，请重新运行脚本，并输入正确物种名称 [human/mouse]";exit
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
     bash $script $read1 10
     reads_length=$(ls *readlengthcount*)
     python3 ~/Tanyuanqing/scripts/XR-polyA-seq_scripts/readslength_distribution_fig_for_XRseq.py $reads_length
     bed_seq=$(ls *rmALL*)
     python3 ~/Tanyuanqing/scripts/XR-polyA-seq_scripts/XRseq_reads_2nt_counts.py $bed_seq
     count_2nt=$(ls *.csv)
     python3 ~/Tanyuanqing/scripts/XR-polyA-seq_scripts/XRseq_reads_2nt_counts_plot.py $count_2nt 
done
