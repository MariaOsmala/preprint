# PRobabilistic Enhancer PRedictIoN Tool PREPRINT: Data preprocessing steps

## Align the reads

```
path_to_dir=.../preprint

data_folder=$path_to_dir/Data

BOWTIE=$path_to_dir/softwares/bowtie2-2.3.3.1/bowtie2
BOWTIE2_INDEXES=$path_to_dir/softwares/genome_indexes/Homo_sapiens/UCSC/hg19/Sequence/Bowtie2Index/genome
BOWTIE_PARAMETERS="-x "$BOWTIE2_INDEXES
cell_line=K562
RAW_data=$path_to_dir/Data/$cell_line/raw_data/

BAM_output=$path_to_dir/Data/$cell_line/bam_replicates/


cd $RAW_data
file=wgEncodeBroadHistoneK562CtcfStdRawDataRep1.fastq.gz
$BOWTIE $BOWTIE_PARAMETERS -U $file | samtools view -bS - | samtools sort - | samtools rmdup -s - -  > $BAM_output${file%".fastq"*}".bam" 
samtools index $BAM_output${file%".fastq"*}".bam"
```

## Process the bam files

Combine bam files of replicates, sort and make indexes. Example for H2az ChIP-seq data in cell line K562.

```
cell_line=K562
data_type=H2az

output_folder=$path_to_dir/Data/$cell_line/bam_combined/

#combine bam files
cd $path_to_dir/Data/$cell_line/bam_replicates/

bam_files=($(ls *${data_type}*.bam))

number_of_files=1

if [ "${#bam_files[@]}" -eq "${number_of_files}" ]; then
    #only one file, just copy it and its index
    cp $bam_files $output_folder
    bam_index=(*$data_type*.bam.bai)
    cp $bam_index $output_folder

else
    #combine files and create a new index
    samtools merge - $bam_files | samtools sort - > $output_folder""$data_type.bam
    #sort the results and make index
    samtools index $output_folder""$data_type.bam
    
fi

	
```

Rename some files: remove cell line names and Rep from files which do not have replicates, so that the names contain only the modification.

## Estimate the fragment lengths using Phantompeakqualtools, Example for H3K4me1 in cell line K562

```
cell_line=K562
data_type=H3k4me1
input_file=Control
data_folder=$path_to_dir/Data

run_spp_path=$path_to_dir/softwares/phantompeakqualtools/

cd $path_to_dir/Data/$cell_line"/bam_combined"

bamfile=$data_type".bam"
outFile=$data_type".out"
inputFile=$input_file".bam"
crossCorrFile=$data_type".pdf" 
outDir=$path_to_dir/Data/$cell_line"/phantompeakqualtools/"
RDataFile=$data_type".RData"

R --slave --args -c=$bamfile -p=5 -i=$inputFile -s=-200:1:500 -odir=$outDir -savd=$outDir$RDataFile -savp=$outDir$crossCorrFile -out=$outDir$outFile -rf < $run_spp_path""run_spp.R

```

## Shift the reads 

```phantompeakqualtools.txt``` contains the estimated shifts and is generated with 
```
R --slave --args --pathToDir=$path_to_dir < code/quality_control_summary.R 
```
Input and DNase-seq are not shifted. For MNase-seq data the shift is 149. This shifts are divided by 2 and rounded to the nearest integer.

For shifting, bam format is converted to bed format, example for Ctcf in cell line K562
```
cell_line=K562
data_type=Ctcf
shift=183
data_folder=$path_to_dir"/Data/"

#convert bam to bed, shift the reads, convert bed back to bam

bamdir=$path_to_dir/Data/$cell_line"/bam_combined/"
bamfile=$data_type".bam"

beddir=$path_to_dir/Data/$cell_line"/bed_combined/"
bedfile=${bamfile%"."*}".bed"


cd $bamdir
bedtools bamtobed -i $bamfile > $beddir$bedfile

genome_folder=$path_to_dir/softwares/bedtools2/genomes/

cd $beddir

beddir_shifted=$path_to_dir/Data/$cell_line"/bed_shifted/"
shift_correct=$( ps -aef | awk  -v shift=$shift 'BEGIN { rounded = sprintf("%.0f", (shift+1)/2); print rounded }' )
shift=shift_correct

shiftBed -i $bedfile -g $genome_folder"human.hg19.genome" -s $shift > $beddir_shifted$bedfile

#convert shift back to bam

cd $beddir_shifted

bamdir_shifted=$path_to_dir/Data/$cell_line"/bam_shifted/"

bedtools bedtobam -i $bedfile -g $genome_folder"human.hg19.genome" | samtools sort - > $bamdir_shifted$bamfile

samtools index $bamdir_shifted$bamfile
``` 


## Convert bed files to format accepted by RFECS 

The RFECS requires the bed format files to have 5 columns and no chrM data.
Remove whitespaces except tabs and chrM data from bed files (needed for RFECS to work). Example for Control data in K562.

```
cd $path_to_dir/Data/K562/bed_shifted/

for f in `ls`; do

    echo $f
    sed 's/ \+//g' $f > "../bed_shifted_RFECS/"$f
    awk -F '\t'  'BEGIN {OFS="\t"} { if (($1) !="chrM")  print }' ../bed_shifted_RFECS/$f > ../bed_shifted_RFECS/tmp
    mv ../bed_shifted_RFECS/tmp "../bed_shifted_RFECS/"$f
    awk -F'\t' 'BEGIN{OFS="\t"}{$5=""; gsub(FS"+",FS); print $0}' "../bed_shifted_RFECS/"$f >  ../bed_shifted_RFECS/tmp
    mv ../bed_shifted_RFECS/tmp "../bed_shifted_RFECS/"$f
done

```
## Next steps

After preprocessing all ChIP-seq, DNase-seq and MNase-seq data for both cell lines, we continue with the [PREPRINT training and genome-wide prediction](TrainingPrediction.md)





