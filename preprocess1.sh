#! /bin/bash

PREPRINT_DIR=`pwd`

data_dir=$PREPRINT_DIR/Data

# Step 1: Align the reads
bowtie_indexes=$PREPRINT_DIR/softwares/genome_indexes/Homo_sapiens/UCSC/hg19/Sequence/Bowtie2Index/genome
bowtie_parameters="-x $bowtie_indexes"

cell_line=K562
raw_data_dir=$data_dir/$cell_line/raw_data
bam_output_dir=$data_dir/$cell_line/bam_replicates

raw_file=wgEncodeBroadHistoneK562CtcfStdRawDataRep1.fastq.gz
bam_output_file=${raw_file%".fastq"*}".bam"

bowtie2 $bowtie_parameters -U $raw_data_dir/$raw_file \
	| samtools view -bS - \
	| samtools sort - \
	| samtools rmdup -s - - \
	> $bam_output_dir/$bam_output_file

samtools index $bam_output_dir/$bam_output_file
