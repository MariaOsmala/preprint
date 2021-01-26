#!/bin/bash

chr=$1
cell_line=$2
bin=$3
path_to_dir=$4


path_str=$path_to_dir"/"$cell_line"/intervals_data_"$bin"/"
cd $path_str

#extract bins that are nonzero
awk '{for(i=4; i<=NF;i++) j+=$i; print j; j=0 }' $path_str"all_"$chr".bedGraph" > $path_str"temp_"$chr".txt" #sum of all modifications

sed -i '1s/.*/sum/' $path_str"temp_"$chr".txt" #add "sum" as column name

paste $path_str"all_"$chr".bedGraph" $path_str"temp_"$chr".txt" > $path_str"temp2_"$chr".txt"

awk '$21 == "0" { next } { print }' $path_str"temp2_"$chr".txt" > $path_str"temp_"$chr".txt"

#remove the sum column

awk '!($21="")' $path_str"temp_"$chr".txt" > $path_str"nozero_"$chr".bedGraph"
rm $path_str"temp_"$chr".txt"
rm $path_str"temp2_"$chr".txt"

#regions with nozero signal    
cut -d' ' -f1-3 $path_str"nozero_"$chr".bedGraph" > $path_str"nozero_regions_only_"$chr".bed"
