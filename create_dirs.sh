#!/bin/bash


path_to_dir=../preprint #specify
cd $path_to_dir


mkdir $path_to_dir/softwares
mkdir $path_to_dir/softwares/R_packages
mkdir $path_to_dir/softwares/R_packages/R-3.4.3
mkdir $path_to_dir/softwares/utils
mkdir $path_to_dir/softwares/bin
mkdir -p $path_to_dir/softwares/genome_indexes/Homo_sapiens/UCSC/hg19/Sequence/Bowtie2Index/genome # put BOWTIE2 indexes here

mkdir $path_to_dir/Data
mkdir $path_to_dir/Data/K562
mkdir $path_to_dir/Data/GM12878
mkdir $path_to_dir/Data/blacklists
mkdir $path_to_dir/Data/GENCODE_TSS
mkdir $path_to_dir/Data/K562/ChromHMM
mkdir $path_to_dir/Data/GM12878/ChromHMM
mkdir $path_to_dir/Data/K562/TFs #narrowPeak files
mkdir $path_to_dir/Data/GM12878/TFs #narrowPeak files
mkdir $path_to_dir/Data/K562/raw_data #fastq files
mkdir $path_to_dir/Data/GM12878/raw_data #fastq files
mkdir $path_to_dir/Data/K562/bam_replicates #DNase-seq and MNase-seq data as .bam and .bam.bai here
mkdir $path_to_dir/Data/GM12878/bam_replicates #DNase-seq and MNase-seq data as .bam and .bam.bai here

#put p300 peak files to these folders
#$path_to_dir/Data/K562/TFs/wgEncodeAwgTfbsSydhK562P300IggrabUniPk.narrowPeak.gz
#$path_to_dir/Data/GM12878/TFs/wgEncodeAwgTfbsSydhGm12878P300IggmusUniPk.narrowPeak.gz

#and DNase-seq peak files to these folders
#$path_to_dir/Data/K562/raw_data/wgEncodeOpenChromDnaseK562PkV2.narrowPeak.gz
#$path_to_dir/Data/GM12878/raw_data/wgEncodeOpenChromDnaseGm12878Pk.narrowPeak.gz

mkdir $path_to_dir/Data/K562/bam_combined
mkdir $path_to_dir/Data/GM12878/bam_combined
mkdir $path_to_dir/Data/K562/phantompeakqualtools
mkdir $path_to_dir/Data/GM12878/phantompeakqualtools

mkdir $path_to_dir/Data/K562/bed_combined
mkdir $path_to_dir/Data/K562/bed_shifted
mkdir $path_to_dir/Data/K562/bam_shifted

mkdir $path_to_dir/Data/GM12878/bed_combined
mkdir $path_to_dir/Data/GM12878/bed_shifted
mkdir $path_to_dir/Data/GM12878/bam_shifted

mkdir $path_to_dir/Data/K562/bed_shifted_RFECS
mkdir $path_to_dir/Data/GM12878/bed_shifted_RFECS

mkdir $path_to_dir/Data/K562/coverage_bigWig
mkdir $path_to_dir/Data/GM12878/coverage_bigWig


mkdir $path_to_dir/figures
mkdir $path_to_dir/Data/K562/data_R
mkdir $path_to_dir/Data/GM12878/data_R

mkdir $path_to_dir/Data/intervals_bed_100
#intervals_data_folder, data in 100 bp intervals, whole genome
mkdir $path_to_dir/Data/K562/intervals_data_100
mkdir $path_to_dir/Data/GM12878/intervals_data_100


mkdir $path_to_dir/results
mkdir $path_to_dir/results/model_promoters_and_random_combined
mkdir $path_to_dir/results/model_promoters_and_random_combined/K562
mkdir $path_to_dir/results/model_promoters_and_random_combined/GM12878

mkdir $path_to_dir/results/model_promoters_and_random_combined/K562/ML
mkdir $path_to_dir/results/model_promoters_and_random_combined/K562/ML/5-fold_CV_1
mkdir $path_to_dir/results/model_promoters_and_random_combined/K562/ML/5-fold_CV_2
mkdir $path_to_dir/results/model_promoters_and_random_combined/K562/ML/5-fold_CV_3
mkdir $path_to_dir/results/model_promoters_and_random_combined/K562/ML/5-fold_CV_4
mkdir $path_to_dir/results/model_promoters_and_random_combined/K562/ML/5-fold_CV_5
mkdir $path_to_dir/results/model_promoters_and_random_combined/K562/Bayes_estimated_priors
mkdir $path_to_dir/results/model_promoters_and_random_combined/K562/Bayes_estimated_priors/5-fold_CV_1
mkdir $path_to_dir/results/model_promoters_and_random_combined/K562/Bayes_estimated_priors/5-fold_CV_2
mkdir $path_to_dir/results/model_promoters_and_random_combined/K562/Bayes_estimated_priors/5-fold_CV_3
mkdir $path_to_dir/results/model_promoters_and_random_combined/K562/Bayes_estimated_priors/5-fold_CV_4
mkdir $path_to_dir/results/model_promoters_and_random_combined/K562/Bayes_estimated_priors/5-fold_CV_5

mkdir $path_to_dir/results/model_promoters_and_random_combined/GM12878/ML
mkdir $path_to_dir/results/model_promoters_and_random_combined/GM12878/Bayes_estimated_priors

mkdir $path_to_dir/results/RFECS_combined
mkdir $path_to_dir/results/RFECS_combined/K562
mkdir $path_to_dir/results/RFECS_combined/K562/training
mkdir $path_to_dir/results/RFECS_combined/K562/trained_forests
mkdir $path_to_dir/results/RFECS_combined/K562/whole_genome_predictions
mkdir $path_to_dir/results/RFECS_combined/K562/bedfiles

mkdir $path_to_dir/results/RFECS_combined/GM12878

mkdir $path_to_dir/results/RFECS_combined/GM12878/whole_genome_predictions
mkdir $path_to_dir/results/RFECS_combined/GM12878/bedfiles



