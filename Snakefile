from snakemake.exceptions import WorkflowError

# The cell line to process
cell_line = 'K562'

# Various paths we will be using in this analysis pipeline
data_dir = 'Data'
softwares_dir = 'softwares'
raw_data_dir = f'{data_dir}/{cell_line}/raw_data'
bam_replicates_dir = f'{data_dir}/{cell_line}/bam_replicates'
bam_combined_dir = f'{data_dir}/{cell_line}/bam_combined'
bam_phantompeakqualtools_dir = f'{data_dir}/{cell_line}/phantompeakqualtools'
bam_shifted_dir = f'{data_dir}/{cell_line}/bam_shifted'
bed_combined_dir = f'{data_dir}/{cell_line}/bed_combined'
bed_shifted_dir = f'{data_dir}/{cell_line}/bed_shifted'

# Used for shifting with bedtools
genome_file = 'bedtools_genomes/human.hg19.genome'

# Step 1: Align the reads
bowtie_indexes = f'{softwares_dir}/genome_indexes/Homo_sapiens/UCSC/hg19/Sequence/Bowtie2Index/genome'
rule align_reads:
	input:
		f'{raw_data_dir}/{{sample}}.fastq.gz'
	output:
		f'{bam_replicates_dir}/{{sample}}.bam'
	shell:
		r'''
		bowtie2 -x {bowtie_indexes} -U {input} \
		| samtools view -bS - \
		| samtools sort - \
		| samtools rmdup -s - - \
		> {output}
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
		bam = f'{bam_combined_dir}/{{data_type}}.bam',
		bai = f'{bam_combined_dir}/{{data_type}}.bam.bai'
	run:
		if len(input) == 0:
			raise WorkflowError(f'No raw files found for data type "{data_type}".')
		elif len(input) == 2:
			# Just copy over the input files
			shell('cp {input[0]} {output.bam}')
			shell('cp {input[1]} {output.bai}')
		else:
			# Merge the input .bam files
			input_bam_files = ' '.join(input[::2])
			shell(
				f'''
				samtools merge - {input_bam_files} \
				| samtools sort - \
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
		r'''
		run_spp.R \
			-c={input.bam} \
			-p=5 \
			-i={input.control} \
			-s=-200:1:500 \
			-rf \
			-odir={bam_phantompeakqualtools_dir} \
			-savd={output.savd} \
			-savp={output.savp} \
			-out={output.out}
		'''

# Step 4: Shift the reads
rule estimate_shifts:
	input:
		'code/quality_control_summary.R'
	output:
		f'{data_dir}/phantompeakqualtools.txt'
	shell:
		'Rscript code/quality_control_summary.R --datadir={data_dir}'

# For shifting, bam format is converted to bed format
rule bam_to_bed:
	input: f'{bam_combined_dir}/{{data_type}}.bam'
	output: f'{bed_combined_dir}/{{data_type}}.bed'
	shell: 'bedtools bamtobed -i {input} > {output}'

rule bed_to_bam:
	input: f'{bed_shifted_dir}/{{data_type}}.bed'
	output: f'{bam_shifted_dir}/{{data_type}}.bam'
	shell:
		r'''
		bedtools bedtobam -i {input} -g {genome_file} \
		| samtools sort - \
		> {output}
		'''

# Input and DNase-seq are not shifted. For MNase-seq data the shift is 149. This shifts are divided by 2 and rounded to the nearest integer.
shift = round(183 / 2)
rule shift_reads:
	input: f'{bed_combined_dir}/{{data_type}}.bed'
	output: f'{bed_shifted_dir}/{{data_type}}.bed'
	run:
		import pandas as pd
		shifts = pd.read_csv(f'{data_dir}/phantompeakqualtools.txt', sep='\t', index_col=[0, 1], usecols=[0, 1, 3])
		shift = shifts.loc[(f'{cell_line}', f'{wildcards.data_type}.bam')][0]
		shift = round(shift / 2)
		shell(f'shiftBed -i {input} -g {genome_file} -s {shift} > {output}')
