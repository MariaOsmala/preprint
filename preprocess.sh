#! /bin/bash

PREPRINT_DIR=`pwd`
DATA_DIR=$PREPRINT_DIR/Data

# Step 1: Align the reads
BOWTIE_INDEXES=$PREPRINT_DIR/softwares/genome_indexes/Homo_sapiens/UCSC/hg19/Sequence/Bowtie2Index/genome
BOWTIE_PARAMETERS="-x $BOWTIE_INDEXES"

CELL_LINE=K562
RAW_DATA_DIR=$PREPRINT_DIR/Data/$CELL_LINE/raw_data
BAM_OUTPUT_DIR=$PREPRINT_DIR/Data/$CELL_LINE/bam_replicates

RAW_FILE=$RAW_DATA_DIR/wgEncodeBroadHistoneK562CtcfStdRawDataRep1.fastq.gz
BAM_OUTPUT_FILE=$BAM_OUTPUT/${RAW_FILE%".fastq"*}".bam"

echo bowtie2 $BOWTIE_PARAMETERS -U $RAW_FILE

bowtie $BOWTIE_PARAMETERS -U $RAW_FILE \
	| samtools view -bS - \
	| samtools sort - \
	| samtools rmdup -s - - \
	> $BAM_OUTPUT_FILE

samtools index $BAM_OUTPUT_FILE
