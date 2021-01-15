# The cell line to process
cell_line = 'K562'

# Various paths we will be using in this analysis pipeline
data_dir = 'Data'
softwares_dir = 'softwares'
raw_data_dir = f'{data_dir}/{cell_line}/raw_data'
bam_replicates_dir = f'{data_dir}/{cell_line}/bam_replicates'
bam_combined_dir = f'{data_dir}/{cell_line}/bam_combined'
bam_phantompeakqualtools_dir = f'{data_dir}/{cell_line}/phantompeakqualtools'

# Step 1: Align the reads
bowtie_indexes = f'{softwares_dir}/genome_indexes/Homo_sapiens/UCSC/hg19/Sequence/Bowtie2Index/genome'
rule align_reads:
	input:
		f'{raw_data_dir}/{{sample}}.fastq.gz'
	output:
		f'{bam_replicates_dir}/{{sample}}.bam'
	shell:
		f'''
		bowtie2 -x {bowtie_indexes} -U {{input}} \\
		| samtools view -bS - \\
		| samtools sort - \\
		| samtools rmdup -s - - \\
		> {{output}}
		'''

# Generic rule to make indexes
rule make_index:
	input: '{sample}.bam'
	output: '{sample}.bam.bai'
	shell: 'samtools index {input}'


# Step 2: Combine bam files of replicates, sort and make indexes
def find_replicates(wildcards):
	"""
	To figure out the list of individual .bam files required to produce the
	combined .bam file for a given data type, we first look at which .fastq.gz
	files are present for the data type using the glob pattern
	'{raw_data_dir}/*{data_type}*.fastq.gz', and then translate them to the
	corresponding .bam and .bam.bai filenames.
	"""
	from glob import glob
	from os.path import basename
	fnames = list()
	print('Finding files for', wildcards.data_type)
	for fname in glob(f'{raw_data_dir}/*{wildcards.data_type}*.fastq.gz'):
		sample = basename(fname).replace('.fastq.gz', '')
		fnames.append(f'{bam_replicates_dir}/{sample}.bam')
		fnames.append(f'{bam_replicates_dir}/{sample}.bam.bai')
	print('Found:', fnames)
	return fnames

rule combine_replicates:
	input:
		find_replicates
	output:
		bam = f'{bam_combined_dir}/{{data_type}}.bam',
		bai = f'{bam_combined_dir}/{{data_type}}.bam.bai'
	run:
		print('Inputs:', input)
		if len(input) == 2:
			# Just copy over the input
			shell('cp {input[0]} {output.bam}')
			shell('cp {input[1]} {output.bai}')
		else:
			input_bam_files = ' '.join(input[::2])
			# Merge the input files
			shell(
				'''
				samtools merge - {input_bam_files} \\
				| samtools sort - \\
				> {output}
				'''
			)
			# Make index for the merged file
			shell('samtools index {output}')


# Step 3: Estimate the fragment lengths using Phantompeakqualtools
rule phantompeakqualtools:
	input:
		bam = f'{bam_combined_dir}/{{data_type}}.bam',
		control = f'{bam_combined_dir}/Control.bam'
	output:
		savd = f'{bam_phantompeakqualtools_dir}/{{data_type}}.RData',
		savp = f'{bam_phantompeakqualtools_dir}/{{data_type}}.pdf',
		out = f'{bam_phantompeakqualtools_dir}/{{data_type}}.out'
	shell:
		f'''
		run_spp.R \\
			-c={{input.bam}} \\
			-p=5 \\
			-i={{input.control}} \\
			-s=-200:1:500 \\
			-rf \\
			-odir={bam_phantompeakqualtools_dir} \\
			-savd={{output.savd}} \\
			-savp={{output.savp}} \\
			-out={{output.out}}
		'''
