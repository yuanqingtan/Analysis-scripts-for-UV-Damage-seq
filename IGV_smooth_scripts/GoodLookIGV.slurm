#!/usr/bin/env bash

#SBATCH -J GoodLookIGV   
#SBATCH -p dna             
#SBATCH -N 1                 
#SBATCH --mem=20G                 
#SBATCH --cpus-per-task=5  
#SBATCH -t 03-00:00:00        
#SBATCH -o GoodLookIGV.Info.out      
#SBATCH --get-user-env

#-------使用说明---------
echo "该脚本用来将单碱基分辨率的reads.bam文件按照自选的binsize和soomthLength生成bigwig文件以获取好看的IGV snapshot图片展示"
echo "用法如下："
echo "sbatch GoodLookIGV.slurm binsize smoothLength"
#-------使用说明---------

binsize=$1
smoothLength=$2

#提取测序数据文件夹名
ls -d Sample* > sn.list

script=~/Tanyuanqing/scripts/IGV_GoodLook/DamageBam_to_bigwig_with_designatedBinsize_byDeeptools_for_goodlookIGV.sh

#获取绝对路径
path=$PWD

#遍历每个数据文件夹里的read1和read2并进行脚本分析
for sl in $(cat sn.list)
do
     cd $path/$sl/
     bam=$(ls *.bam)
     bash $script $bam $binsize $smoothLength
done
