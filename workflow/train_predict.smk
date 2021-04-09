# Rules for downloading various files
rule download_gencode:
	input:
	output: f'{gencode_dir}/gencode.v27lift37.annotation.gtf.gz'
	shell: 'wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_27/GRCh37_mapping/gencode.v27lift37.annotation.gtf.gz -O {output}'

rule download_blacklists:
	input:
	output: f'{blacklists_dir}/{{file}}'
	shell: 'wget http://hgdownload.cse.ucsc.edu/goldenpath/hg19/encodeDCC/wgEncodeMapability/{wildcards.file} -O {output}'

# Function that finds all the available BAM files
def all_bam_files(wildcards):
	return [
		f'{data_dir}/{wildcards.cell_line}/bam_shifted/{data_type}.bam'
		for data_type in all_data_types(wildcards.cell_line)
	]

# The steps needed for the analysis
rule whole_genome_coverage:
	input:
		code=f'{scripts_dir}/1_whole_genome_coverage.R',
		bam_files=all_bam_files,
	output:
		f'{data_dir}/{{cell_line}}/data_R/whole_genome_coverage.rds'
	shell:
		f'Rscript {scripts_dir}/1_whole_genome_coverage.R'

rule make_profiles:
	input:
		code=f'{scripts_dir}/2_make_profiles.R',
		bam_files=all_bam_files,
		p300=f'{data_dir}/{{cell_line}}/raw_data/wgEncodeAwgTfbsSydhK562P300IggrabUniPk.narrowPeak.gz',
		DNase=f'{data_dir}/{{cell_line}}/raw_data/wgEncodeOpenChromDnaseK562PkV2.narrowPeak.gz',
		blacklist_Dac=f'{blacklists_dir}/wgEncodeDacMapabilityConsensusExcludable.bed.gz',
		blacklist_Duke=f'{blacklists_dir}/wgEncodeDukeMapabilityRegionsExcludable.bed.gz',
		TSS_annotation=f'{gencode_dir}/gencode.v27lift37.annotation.gtf.gz',
		whole_genome_cov=f'{data_dir}/{{cell_line}}/data_R/whole_genome_coverage.rds',
	output:
		f'{data_dir}/{{cell_line}}/data_R/profiles.rds'
	shell:
		f'Rscript {scripts_dir}/2_make_profiles.R'

rule train_predict:
	input:
		code=f'{scripts_dir}/3_train_predict.R',
		profiles=f'{data_dir}/{{cell_line}}/data_R/profiles.rds',
	output:
		f'{data_dir}/{{cell_line}}/data_R/predictions.RData'
	shell:
		f'Rscript {scripts_dir}/3_train_predict.R'
