# The cell line to process
cell_line = 'K562'

# Various paths we will be using in this analysis pipeline
data_dir = 'Data'
softwares_dir = 'softwares'
raw_data_dir = f'{data_dir}/{cell_line}/raw_data'
bam_replicates_dir = f'{data_dir}/{cell_line}/bam_replicates'
bam_combined_dir = f'{data_dir}/{cell_line}/bam_combined'

print(raw_data_dir)
print(bam_replicates_dir)

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
	for fname in glob(f'{raw_data_dir}/*{wildcards.data_type}*.fastq.gz'):
		sample = basename(fname).replace('.fastq.gz', '')
		fnames.append(f'{bam_replicates_dir}/{sample}.bam')
		fnames.append(f'{bam_replicates_dir}/{sample}.bam.bai')
	return fnames

rule combine_replicates:
	input:
		find_replicates
	output:
		f'{bam_combined_dir}/{{data_type}}.bam'
	run:
		print(input)
		if len(input) == 1:
			# Just copy over the input
			shell('cp {input} {output}')
			shell('cp {input}.bai {output}.bai')
		else:
			# Merge the input files
			shell(
				'''
				samtools merge - {input} \\
				| samtools sort - \\
				> {output}
				'''
			)
			shell('samtools index {output}')
