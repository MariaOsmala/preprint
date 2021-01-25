rule define_TSS:
	input:
		f'{gencode_dir}/gencode.v27lift37.annotation.gtf.gz',
		f'{code_dir}/define_TSS.R'
	output:
		f'{gencode_dir}/GENCODE.RData',
		f'{gencode_dir}/GR_Gencode_protein_coding_TSS.RDS',
		f'{gencode_dir}/GR_Gencode_protein_coding_TSS_positive.RDS',
		f'{gencode_dir}/GR_Gencode_TSS.RDS',
		f'{gencode_dir}/GR_Gencode_TSS_positive.RDS'
	shell:
		'Rscript {code_dir}/define_TSS.R --pathToDir={data_dir}'

rule download_gencode:
	input:
	output: f'{gencode_dir}/gencode.v27lift37.annotation.gtf.gz'
	shell: 'wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_27/GRCh37_mapping/gencode.v27lift37.annotation.gtf.gz -O {output}'

rule download_blacklists:
	input:
	output: f'{blacklists_dir}/{{file}}'
	shell: 'wget http://hgdownload.cse.ucsc.edu/goldenpath/hg19/encodeDCC/wgEncodeMapability/{wildcards.file} -O {output}'

rule extract_enhancers:
	input:
		p300=f'{raw_data_dir}/wgEncodeAwgTfbsSydhK562P300IggrabUniPk.narrowPeak.gz',
		DNase=f'{raw_data_dir}/wgEncodeOpenChromDnaseK562PkV2.narrowPeak.gz',
		blacklist_Dac=f'{blacklists_dir}/wgEncodeDacMapabilityConsensusExcludable.bed.gz',
		blacklist_Duke=f'{blacklists_dir}/wgEncodeDukeMapabilityRegionsExcludable.bed.gz'
	output:
		f'{data_r_dir}/{config["extract_enhancers"]["N"]}_enhancers_bin_{config["extract_enhancers"]["binSize"]}_window_{config["extract_enhancers"]["window"]}.RData'
	shell:
		r'''
		Rscript code/extract_enhancers.R \
			--window={config[extract_enhancers][window]} \
			--binSize={config[extract_enhancers][binSize]} \
			--N={config[extract_enhancers][N]} \
			--distToPromoter={config[extract_enhancers][distToPromotor]} \
			--pathToDir={data_dir} \
			--cellLine={cell_line} \
			--p300File={input.p300} \
			--DNaseFile={input.DNase} \
			--normalize=FALSE \
			--NormCellLine=""
		'''
