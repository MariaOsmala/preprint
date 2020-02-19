#!/bin/bash

#path_to_dir=/m/cs/scratch/csb/projects/enhancer_prediction_BMCBioinformatics/preprint #specify
path_to_dir=/m/cs/scratch/csb/projects/enhancer_prediction/experiments/RProjects/preprint
cd $path_to_dir


mkdir $path_to_dir/softwares/
mkdir $path_to_dir/softwares/R_packages/
mkdir $path_to_dir/softwares/R_packages/R-3.4.3/
mkdir $path_to_dir/softwares/utils/
mkdir $path_to_dir/softwares/bin/
mkdir $path_to_dir/softwares/genome_indexes/.../genome # put BOWTIE2 indexes here

mkdir $path_to_dir/Data
cd $path_to_dir/Data
mkdir K562
mkdir GM12878
mkdir blacklists
mkdir GENCODE_TSS
mkdir K562/ChromHMM
mkdir GM12878/ChromHMM
mkdir K562/TFs #narrowPeak files
mkdir GM12878/TFs #narrowPeak files
mkdir K562/raw_data #fastq files
mkdir GM12878/raw_data #fastq files
mkdir K562/bam_replicates #DNase-seq and MNase-seq data as .bam and .bam.bai here
mkdir GM12878/bam_replicates #DNase-seq and MNase-seq data as .bam and .bam.bai here

#put p300 peak files to these folders
#$path_to_dir/Data/K562/TFs/wgEncodeAwgTfbsSydhK562P300IggrabUniPk.narrowPeak.gz
#$path_to_dir/Data/GM12878/TFs/wgEncodeAwgTfbsSydhGm12878P300IggmusUniPk.narrowPeak.gz

#and DNase-seq peak files to these folders
#$path_to_dir/Data/K562/raw_data/wgEncodeOpenChromDnaseK562PkV2.narrowPeak.gz
#$path_to_dir/Data/GM12878/raw_data/wgEncodeOpenChromDnaseGm12878Pk.narrowPeak.gz

mkdir K562/bam_combined
mkdir GM12878/bam_combined
mkdir K562/phantompeakqualtools
mkdir GM12878/phantompeakqualtools

mkdir $path_to_dir/Data/K562/bed_combined/
mkdir $path_to_dir/Data/K562/bed_shifted/
mkdir $path_to_dir/Data/K562/bam_shifted/

mkdir $path_to_dir/Data/GM12878/bed_combined/
mkdir $path_to_dir/Data/GM12878/bed_shifted/
mkdir $path_to_dir/Data/GM12878/bam_shifted/

mkdir $path_to_dir/Data/K562/bed_shifted_RFECS/
mkdir $path_to_dir/Data/GM12878/bed_shifted_RFECS/

mkdir $path_to_dir/Data/K562/coverage_bigWig/
mkdir $path_to_dir/Data/GM12878/coverage_bigWig/


mkdir $path_to_dir/figures/
mkdir $path_to_dir/Data/K562/data_R/
mkdir $path_to_dir/Data/GM12878/data_R/

mkdir $path_to_dir/results/
mkdir $path_to_dir/results/model_promoters_and_random/
mkdir $path_to_dir/results/do_not_model_promoters_and_random/

mkdir $path_to_dir/results/model_promoters_and_random/K562/
mkdir $path_to_dir/results/model_promoters_and_random/GM12878/
mkdir $path_to_dir/results/do_not_model_promoters_and_random/K562/
mkdir $path_to_dir/results/do_not_model_promoters_and_random/GM12878/

mkdir $path_to_dir/results/model_promoters_and_random/K562/pure_random/
mkdir $path_to_dir/results/model_promoters_and_random/K562/pure_random/ML
mkdir $path_to_dir/results/model_promoters_and_random/K562/pure_random/ML/5-fold_CV_1
mkdir $path_to_dir/results/model_promoters_and_random/K562/pure_random/ML/5-fold_CV_2
mkdir $path_to_dir/results/model_promoters_and_random/K562/pure_random/ML/5-fold_CV_3
mkdir $path_to_dir/results/model_promoters_and_random/K562/pure_random/ML/5-fold_CV_4
mkdir $path_to_dir/results/model_promoters_and_random/K562/pure_random/ML/5-fold_CV_5
mkdir $path_to_dir/results/model_promoters_and_random/K562/pure_random/Bayes_estimated_priors
mkdir $path_to_dir/results/model_promoters_and_random/K562/pure_random/Bayes_estimated_priors/5-fold_CV_1
mkdir $path_to_dir/results/model_promoters_and_random/K562/pure_random/Bayes_estimated_priors/5-fold_CV_2
mkdir $path_to_dir/results/model_promoters_and_random/K562/pure_random/Bayes_estimated_priors/5-fold_CV_3
mkdir $path_to_dir/results/model_promoters_and_random/K562/pure_random/Bayes_estimated_priors/5-fold_CV_4
mkdir $path_to_dir/results/model_promoters_and_random/K562/pure_random/Bayes_estimated_priors/5-fold_CV_5
mkdir $path_to_dir/results/model_promoters_and_random/K562/random_with_signal/
mkdir $path_to_dir/results/model_promoters_and_random/K562/random_with_signal/ML
mkdir $path_to_dir/results/model_promoters_and_random/K562/random_with_signal/ML/5-fold_CV_1
mkdir $path_to_dir/results/model_promoters_and_random/K562/random_with_signal/ML/5-fold_CV_2
mkdir $path_to_dir/results/model_promoters_and_random/K562/random_with_signal/ML/5-fold_CV_3
mkdir $path_to_dir/results/model_promoters_and_random/K562/random_with_signal/ML/5-fold_CV_4
mkdir $path_to_dir/results/model_promoters_and_random/K562/random_with_signal/ML/5-fold_CV_5
mkdir $path_to_dir/results/model_promoters_and_random/K562/random_with_signal/Bayes_estimated_priors
mkdir $path_to_dir/results/model_promoters_and_random/K562/random_with_signal/Bayes_estimated_priors/5-fold_CV_1
mkdir $path_to_dir/results/model_promoters_and_random/K562/random_with_signal/Bayes_estimated_priors/5-fold_CV_2
mkdir $path_to_dir/results/model_promoters_and_random/K562/random_with_signal/Bayes_estimated_priors/5-fold_CV_3
mkdir $path_to_dir/results/model_promoters_and_random/K562/random_with_signal/Bayes_estimated_priors/5-fold_CV_4
mkdir $path_to_dir/results/model_promoters_and_random/K562/random_with_signal/Bayes_estimated_priors/5-fold_CV_5
mkdir $path_to_dir/results/model_promoters_and_random/GM12878/pure_random/
mkdir $path_to_dir/results/model_promoters_and_random/GM12878/pure_random/ML
mkdir $path_to_dir/results/model_promoters_and_random/GM12878/pure_random/Bayes_estimated_priors
mkdir $path_to_dir/results/model_promoters_and_random/GM12878/random_with_signal/
mkdir $path_to_dir/results/model_promoters_and_random/GM12878/random_with_signal/ML
mkdir $path_to_dir/results/model_promoters_and_random/GM12878/random_with_signal/Bayes_estimated_priors

mkdir $path_to_dir/results/do_not_model_promoters_and_random/K562/pure_random/
mkdir $path_to_dir/results/do_not_model_promoters_and_random/K562/pure_random/ML
mkdir $path_to_dir/results/do_not_model_promoters_and_random/K562/pure_random/ML/5-fold_CV_1
mkdir $path_to_dir/results/do_not_model_promoters_and_random/K562/pure_random/ML/5-fold_CV_2
mkdir $path_to_dir/results/do_not_model_promoters_and_random/K562/pure_random/ML/5-fold_CV_3
mkdir $path_to_dir/results/do_not_model_promoters_and_random/K562/pure_random/ML/5-fold_CV_4
mkdir $path_to_dir/results/do_not_model_promoters_and_random/K562/pure_random/ML/5-fold_CV_5
mkdir $path_to_dir/results/do_not_model_promoters_and_random/K562/pure_random/Bayes_estimated_priors
mkdir $path_to_dir/results/do_not_model_promoters_and_random/K562/pure_random/Bayes_estimated_priors/5-fold_CV_1
mkdir $path_to_dir/results/do_not_model_promoters_and_random/K562/pure_random/Bayes_estimated_priors/5-fold_CV_2
mkdir $path_to_dir/results/do_not_model_promoters_and_random/K562/pure_random/Bayes_estimated_priors/5-fold_CV_3
mkdir $path_to_dir/results/do_not_model_promoters_and_random/K562/pure_random/Bayes_estimated_priors/5-fold_CV_4
mkdir $path_to_dir/results/do_not_model_promoters_and_random/K562/pure_random/Bayes_estimated_priors/5-fold_CV_5
mkdir $path_to_dir/results/do_not_model_promoters_and_random/K562/random_with_signal/
mkdir $path_to_dir/results/do_not_model_promoters_and_random/K562/random_with_signal/ML
mkdir $path_to_dir/results/do_not_model_promoters_and_random/K562/random_with_signal/ML/5-fold_CV_1
mkdir $path_to_dir/results/do_not_model_promoters_and_random/K562/random_with_signal/ML/5-fold_CV_2
mkdir $path_to_dir/results/do_not_model_promoters_and_random/K562/random_with_signal/ML/5-fold_CV_3
mkdir $path_to_dir/results/do_not_model_promoters_and_random/K562/random_with_signal/ML/5-fold_CV_4
mkdir $path_to_dir/results/do_not_model_promoters_and_random/K562/random_with_signal/ML/5-fold_CV_5
mkdir $path_to_dir/results/do_not_model_promoters_and_random/K562/random_with_signal/Bayes_estimated_priors
mkdir $path_to_dir/results/do_not_model_promoters_and_random/K562/random_with_signal/Bayes_estimated_priors/5-fold_CV_1
mkdir $path_to_dir/results/do_not_model_promoters_and_random/K562/random_with_signal/Bayes_estimated_priors/5-fold_CV_2
mkdir $path_to_dir/results/do_not_model_promoters_and_random/K562/random_with_signal/Bayes_estimated_priors/5-fold_CV_3
mkdir $path_to_dir/results/do_not_model_promoters_and_random/K562/random_with_signal/Bayes_estimated_priors/5-fold_CV_4
mkdir $path_to_dir/results/do_not_model_promoters_and_random/K562/random_with_signal/Bayes_estimated_priors/5-fold_CV_5
mkdir $path_to_dir/results/do_not_model_promoters_and_random/GM12878/pure_random/
mkdir $path_to_dir/results/do_not_model_promoters_and_random/GM12878/pure_random/ML
mkdir $path_to_dir/results/do_not_model_promoters_and_random/GM12878/pure_random/Bayes_estimated_priors
mkdir $path_to_dir/results/do_not_model_promoters_and_random/GM12878/random_with_signal/
mkdir $path_to_dir/results/do_not_model_promoters_and_random/GM12878/random_with_signal/ML
mkdir $path_to_dir/results/do_not_model_promoters_and_random/GM12878/random_with_signal/Bayes_estimated_priors

mkdir $path_to_dir/results/RFECS/
mkdir $path_to_dir/results/RFECS/K562
mkdir $path_to_dir/results/RFECS/K562/pure_random/
mkdir $path_to_dir/results/RFECS/K562/pure_random/training/
mkdir $path_to_dir/results/RFECS/K562/pure_random/trained_forests/
mkdir $path_to_dir/results/RFECS/K562/pure_random/whole_genome_predictions/
mkdir $path_to_dir/results/RFECS/K562/pure_random/bedfiles/
mkdir $path_to_dir/results/RFECS/K562/random_with_signal/
mkdir $path_to_dir/results/RFECS/K562/random_with_signal/training/
mkdir $path_to_dir/results/RFECS/K562/random_with_signal/trained_forests/
mkdir $path_to_dir/results/RFECS/K562/random_with_signal/whole_genome_predictions/
mkdir $path_to_dir/results/RFECS/K562/random_with_signal/bedfiles/
mkdir $path_to_dir/results/RFECS/GM12878
mkdir $path_to_dir/results/RFECS/GM12878/pure_random/

mkdir $path_to_dir/results/RFECS/GM12878/pure_random/whole_genome_predictions
mkdir $path_to_dir/results/RFECS/GM12878/pure_random/bedfiles
mkdir $path_to_dir/results/RFECS/GM12878/random_with_signal/

mkdir $path_to_dir/results/RFECS/GM12878/random_with_signal/whole_genome_predictions
mkdir $path_to_dir/results/RFECS/GM12878/random_with_signal/bedfiles

mkdir $path_to_dir/Data/intervals_bed_100


#intervals_data_folder, data in 100 bp intervals, whole genome
mkdir $path_to_dir/Data/K562/intervals_data_100
mkdir $path_to_dir/Data/GM12878/intervals_data_100






