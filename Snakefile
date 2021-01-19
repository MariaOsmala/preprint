from glob import glob
from os.path import basename
from snakemake.exceptions import WorkflowError
import pandas as pd

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
bed_shifted_RFECS_dir = f'{data_dir}/{cell_line}/bed_shifted_RFECS'

# Big list of all the samples
all_samples = pd.read_csv('samples.tsv', sep='\t').sort_index()
data_types = all_samples.query(f"cell_line=='{cell_line}'").data_type.unique()

# The default rule that will do the entire preprocessing pipeline
rule preprocess:
	input:
		expand(f'{bed_shifted_RFECS_dir}/{{data_type}}.bed', data_type=data_types)


# Step 0: Download the data
rule download:
	input:
	output: f'{raw_data_dir}/{{fname}}'
	message: 'Downloading {wildcards.fname}'
	run:
		URL = all_samples.query(f"cell_line=='{cell_line}' and fname=='{wildcards.fname}'")['URL'].item()
		shell(f'wget {URL} -O {output}')



# Step 1: Align the reads
def find_raw_files(wildcards):
	"""
	Uses the all_samples table to find the raw file(s) corresponding to a given
	sample. These is either a .fastq.gz file or both a .bam and .bam.bai file.
	"""
	fnames = list()
	for url in all_samples.query(f"cell_line=='{cell_line}' and sample=='{wildcards.sample}'")['URL']:
		_, fname = url.rsplit('/', 1)
		fnames.append(f'{raw_data_dir}/{fname}')
	return fnames

def get_ext(fname):
	return fname.rsplit('/', 1)[-1].split('.', 1)[-1]

rule align_reads:
	input:
		find_raw_files
	output:
		bam = f'{bam_replicates_dir}/{{sample}}.bam',
		bai = f'{bam_replicates_dir}/{{sample}}.bam.bai'
	run:
		input_ext = get_ext(input[0])
		if input_ext == 'fastq.gz':
			# .fastq.gz input files need aligning with bowtie
			bowtie_indexes = f'{softwares_dir}/genome_indexes/Homo_sapiens/UCSC/hg19/Sequence/Bowtie2Index/genome'
			shell(
				r'''
				bowtie2 -x {bowtie_indexes} -U {input[0]} \
				| samtools view -bS - \
				| samtools sort - \
				| samtools rmdup -s - - \
				> {output.bam}
				samtools index {output.bam}
				''')
		elif input_ext == 'bam':
			# .bam input files do not need aligning. Just copy them.
			shell('cp {input[0]} {output.bam}')
			shell('cp {input[1]} {output.bai}')
		else:
			raise WorkflowError(f'Unknown file type: "{input[0]}"')


# Step 2: Combine bam files of replicates, sort and make indexes
def find_replicates(wildcards):
	"""
	Uses the all_samples table to find the list of individual .bam files
	required to produce the combined .bam file for a given data type.
	"""
	samples = all_samples.query(f"cell_line=='{cell_line}' and data_type=='{wildcards.data_type}'")['sample']
	fnames = list()
	for sample in samples:
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
			raise WorkflowError(f'No raw files found for data type "{wildcards.data_type}".')
		elif len(input) == 2:
			# Just copy over the input files
			shell('cp {input[0]} {output.bam}')
			shell('cp {input[1]} {output.bai}')
		else:
			# Merge the input .bam files
			input_bam_files = ' '.join(input[::2])
			shell(
				r'''
				samtools merge - {input_bam_files} \
				| samtools sort - \
				> {output.bam}
				''')
			# Make index for the merged file
			shell('samtools index {output.bam}')


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
genome_file = 'bedtools_genomes/human.hg19.genome'

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

rule shift_reads:
	input: f'{bed_combined_dir}/{{data_type}}.bed'
	output: f'{bed_shifted_dir}/{{data_type}}.bed'
	run:
		# Input and DNase-seq are not shifted.
		if 'Input' in input or 'DNase' in input:
			shell(f'cp {input} {output}')
		else:
			# For MNase-seq data the shift is 149.
			if 'MNase' in input:
				shift = 149
			else:
				# Lookup the required shift produced by phantompeakqualtools.
				# These shifts are divided by 2 and rounded to the nearest integer.
				shifts = pd.read_csv(f'{data_dir}/phantompeakqualtools.txt', sep='\t', index_col=[0, 1], usecols=[0, 1, 3])
				shift = shifts.loc[(f'{cell_line}', f'{wildcards.data_type}.bam')][0]
				shift = round(shift / 2)
			shell(f'shiftBed -i {input} -g {genome_file} -s {shift} > {output}')


# Step 5: Convert bed files to format accepted by RFECS
rule convert_bed_for_RFECS:
	input: f'{bed_shifted_dir}/{{data_type}}.bed'
	output: f'{bed_shifted_RFECS_dir}/{{data_type}}.bed'
	shell:
		r"""
		sed 's/ \+//g' {input} > {output}
		awk -F '\t' 'BEGIN {{OFS="\t"}} {{ if (($1) != "chrM")  print }}' {output} > {bed_shifted_RFECS_dir}/tmp
		mv {bed_shifted_RFECS_dir}/tmp {output}
		awk -F'\t' 'BEGIN{{OFS="\t"}}{{$5=""; gsub(FS"+",FS); print $0}}' {output} > {bed_shifted_RFECS_dir}/tmp
		mv {bed_shifted_RFECS_dir}/tmp {output}
		"""
