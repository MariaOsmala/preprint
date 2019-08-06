#!/bin/bash

chr=$1
cell_line=$2
bin=$3
path_to_dir=$4


path_str=$path_to_dir"/Data/"$cell_line"/intervals_data_"$bin"/"
cd $path_str

str1="K562"
str2="GM12878"

#for K562
if [ "$cell_line" == "$str1" ]; then
	bedtools unionbedg -header -i Control/$chr.bed Ctcf/$chr.bed H2az/$chr.bed H3k27ac/$chr.bed H3k27me3/$chr.bed H3k36me3/$chr.bed H3k4me1/$chr.bed H3k4me2/$chr.bed H3k4me3/$chr.bed H3k79me2/$chr.bed H3k9ac/$chr.bed H3k9me3/$chr.bed H4k20me1/$chr.bed InputV2/$chr.bed Nsome/$chr.bed OpenChromDnaseV2/$chr.bed Pol2RawData/$chr.bed -names Control Ctcf H2az H3k27ac H3k27me3 H3k36me3 H3k4me1 H3k4me2 H3k4me3 H3k79me2 H3k9ac H3k9me3 H4k20me1 InputV2 Nsome OpenChromDnaseV2 Pol2RawData > $path_str"all_"$chr".bedGraph"
else
#for GM12878
	bedtools unionbedg -header -i Control/$chr.bed Ctcf/$chr.bed H2az/$chr.bed H3k27ac/$chr.bed H3k27me3/$chr.bed H3k36me3/$chr.bed H3k4me1/$chr.bed H3k4me2/$chr.bed H3k4me3/$chr.bed H3k79me2/$chr.bed H3k9ac/$chr.bed H3k9me3/$chr.bed H4k20me1/$chr.bed Input/$chr.bed Nsome/$chr.bed OpenChromDnase/$chr.bed Pol2RawData/$chr.bed -names Control Ctcf H2az H3k27ac H3k27me3 H3k36me3 H3k4me1 H3k4me2 H3k4me3 H3k79me2 H3k9ac H3k9me3 H4k20me1 Input Nsome OpenChromDnase Pol2RawData > $path_str"all_"$chr".bedGraph"
fi

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
