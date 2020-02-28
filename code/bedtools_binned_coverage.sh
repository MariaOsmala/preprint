#!/bin/bash

chr=$1
mod=$2
bin=$3
cell_line=$4
path_to_dir=$5


echo $chr
echo $mod
echo $bin
echo $cell_line
echo $path_to_dir



intervals_data_folder=$path_to_dir"/Data/"$cell_line"/intervals_data_"$bin
echo $intervals_data_folder"/"$mod
mkdir $intervals_data_folder"/"$mod
intervals=$path_to_dir"/Data/intervals_bed_"$bin"/"$chr".bed" 

bam_folder=$path_to_dir"/Data/"$cell_line"/bam_shifted/"
bam_file=$bam_folder""$mod".bam"

#genomeCoverageBed
cd $bam_folder
#By default, overlaps are reported without respect to strand.

bedtools multicov -bams $bam_file -bed $intervals | sort -k 1,1 -k2,2 -n | cut -f 1-3,7 > $intervals_data_folder"/"$mod"/"$chr".bed"

