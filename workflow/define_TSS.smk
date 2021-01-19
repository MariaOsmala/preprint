rule define_TSS:
	input:
		f'{gencode_dir}/gencode.v27lift37.annotation.gtf.gz',
		'code/define_TSS.R'
	output:
		f'{gencode_dir}/GENCODE.RData',
		f'{gencode_dir}/GR_Gencode_protein_coding_TSS.RDS',
		f'{gencode_dir}/GR_Gencode_protein_coding_TSS_positive.RDS',
		f'{gencode_dir}/GR_Gencode_TSS.RDS',
		f'{gencode_dir}/GR_Gencode_TSS_positive.RDS'
	shell:
		'Rscript code/define_TSS.R --pathToDir={data_dir}'

rule download_gencode:
	input:
	output: f'{gencode_dir}/gencode.v27lift37.annotation.gtf.gz'
	shell: 'wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_27/GRCh37_mapping/gencode.v27lift37.annotation.gtf.gz -O {output}'
