from glob import glob
from os.path import basename
from snakemake.exceptions import WorkflowError

# The default rule that will do the entire preprocessing pipeline
rule preprocess:
	input:
		expand(f'{data_dir}/K562/bam_shifted/{{data_type}}.bam', data_type=all_data_types('K562')),
		expand(f'{data_dir}/Gm12878/bam_shifted/{{data_type}}.bam', data_type=all_data_types('Gm12878')),


# Step 0: Download the data
rule download:
	input:
	output: f'{data_dir}/{{cell_line}}/raw_data/{{fname, [^/]+}}'
	message: 'Downloading {wildcards.fname}'
	run:
		URL = all_samples.query(f"cell_line=='{wildcards.cell_line}' and fname=='{wildcards.fname}'")['URL'].item()
		shell(f'wget {URL} -O {output}')


rule download_bowtie2_index:
	input:
	output: f'{data_dir}/bowtie2-index/hg19.zip'
	shell: 'wget https://genome-idx.s3.amazonaws.com/bt/hg19.zip -O {output}'


rule unzip_bowtie2_index:
	input:
	    file = f'{data_dir}/bowtie2-index/hg19.zip'
	output:
	    folder = directory(f'{data_dir}/bowtie2-index/hg19/'),
		output1 = f'{data_dir}/bowtie2-index/hg19/hg19.1.bt2',
		output2 = f'{data_dir}/bowtie2-index/hg19/hg19.2.bt2',
		output3 = f'{data_dir}/bowtie2-index/hg19/hg19.3.bt2',
		output4 = f'{data_dir}/bowtie2-index/hg19/hg19.4.bt2',
		output5 = f'{data_dir}/bowtie2-index/hg19/hg19.rev.1.bt2',
		output6 = f'{data_dir}/bowtie2-index/hg19/hg19.rev.2.bt2',
		output7 = f'{data_dir}/bowtie2-index/hg19/make_hg19.sh',
	shell: 'unzip {input} -d {output.folder} -x {output.output1} {output.output2} {output.output3} {output.output4} {output.output5} {output.output6} {output.output7}'





# Step 1: Align the reads
def find_raw_files(wildcards):
	"""
	Uses the all_samples table to find the raw file(s) corresponding to a given
	sample. These is either a .fastq.gz file or both a .bam and .bam.bai file.
	"""
	fnames = list()
	for url in all_samples.query(f"cell_line=='{wildcards.cell_line}' and sample=='{wildcards.sample}'")['URL']:
		_, fname = url.rsplit('/', 1)
		fnames.append(f'{data_dir}/{wildcards.cell_line}/raw_data/{fname}')
	return fnames

def get_ext(fname):
	return fname.rsplit('/', 1)[-1].split('.', 1)[-1]

rule align_reads:
	input:
		find_raw_files,
		index1 = f'{data_dir}/bowtie2-index/hg19/hg19.1.bt2',
		index2 = f'{data_dir}/bowtie2-index/hg19/hg19.2.bt2',
		index3 = f'{data_dir}/bowtie2-index/hg19/hg19.3.bt2',
		index4 = f'{data_dir}/bowtie2-index/hg19/hg19.4.bt2',
		index5 = f'{data_dir}/bowtie2-index/hg19/hg19.rev.1.bt2',
		index6 = f'{data_dir}/bowtie2-index/hg19/hg19.rev.2.bt2',
		index7 = f'{data_dir}/bowtie2-index/hg19/make_hg19.sh'
	output:
		bam = f'{data_dir}/{{cell_line}}/bam_replicates/{{sample}}.bam',
		bai = f'{data_dir}/{{cell_line}}/bam_replicates/{{sample}}.bam.bai',
   threads: 2
	run:
		input_ext = get_ext(input[0])
		if input_ext == 'fastq.gz':
			# .fastq.gz input files need aligning with bowtie
			bowtie_indexes = f'{data_dir}/bowtie2-index/hg19/hg19'
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
	samples = all_samples.query(f"cell_line=='{wildcards.cell_line}' and data_type=='{wildcards.data_type}'")['sample']
	fnames = list()
	for sample in samples:
		fnames.append(f'{data_dir}/{wildcards.cell_line}/bam_replicates/{sample}.bam')
		fnames.append(f'{data_dir}/{wildcards.cell_line}/bam_replicates/{sample}.bam.bai')
	return fnames

rule combine_replicates:
	input:
		find_replicates
	output:
		bam = f'{data_dir}/{{cell_line}}/bam_combined/{{data_type}}.bam',
		bai = f'{data_dir}/{{cell_line}}/bam_combined/{{data_type}}.bam.bai'
	threads: 5
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
			step_1_threads = threads-1
			step_2_threads = 1
			shell(
				r'''
				samtools merge - {input_bam_files} --threads {step_1_threads} \
				| samtools sort -T /tmp -@ {step_2_threads} - \
				> {output.bam}
				''')
			# Make index for the merged file
			shell('samtools index {output.bam}')


# Step 3: Estimate the fragment lengths using Phantompeakqualtools
rule phantompeakqualtools:
	input:
		bam = f'{data_dir}/{{cell_line}}/bam_combined/{{data_type}}.bam',
		control = f'{data_dir}/{{cell_line}}/bam_combined/Control.bam'
	output:
		savd = f'{data_dir}/{{cell_line}}/phantompeakqualtools/{{data_type}}.RData',
		savp = f'{data_dir}/{{cell_line}}/phantompeakqualtools/{{data_type}}.pdf',
		out = f'{data_dir}/{{cell_line}}/phantompeakqualtools/{{data_type}}.out'
	threads: 4
	shell:
		r'''
		run_spp.R \
			-c={input.bam} \
			-p={threads} \
			-i={input.control} \
			-s=-200:1:500 \
			-rf \
			-odir={data_dir}/{wildcards.cell_line}/phantompeakqualtools \
			-savd={output.savd} \
			-savp={output.savp} \
			-out={output.out}
		'''

# Step 4: Shift the reads
genome_file = f'{project_dir}/bedtools_genomes/human.hg19.genome'

rule estimate_shifts:
	input:
		f'{scripts_dir}/quality_control_summary.R',
		expand(f'{data_dir}/K562/phantompeakqualtools/{{data_type}}.out', data_type=all_data_types('K562', exclude=['OpenChromDnaseV2', 'InputV2', 'Nsome'])),
		expand(f'{data_dir}/Gm12878/phantompeakqualtools/{{data_type}}.out', data_type=all_data_types('Gm12878', exclude=['OpenChromDnase', 'Input', 'Nsome'])),
	output:
		f'{data_dir}/phantompeakqualtools.txt'
	shell:
		'Rscript {scripts_dir}/quality_control_summary.R --datadir={data_dir}'

# For shifting, bam format is converted to bed format
rule bam_to_bed:
	input: f'{data_dir}/{{cell_line}}/bam_combined/{{data_type}}.bam'
	output: f'{data_dir}/{{cell_line}}/bed_combined/{{data_type}}.bed'
	shell: 'bedtools bamtobed -i {input} > {output}'

rule bed_to_bam:
	input:
		f'{data_dir}/{{cell_line}}/bed_shifted/{{data_type}}.bed'
	output:
		bam=f'{data_dir}/{{cell_line}}/bam_shifted/{{data_type}}.bam',
		bai=f'{data_dir}/{{cell_line}}/bam_shifted/{{data_type}}.bam.bai',
	threads: 2
	shell:
		r'''
		bedtools bedtobam -i {input} -g {genome_file} \
		| samtools sort - \
		> {output.bam}
		samtools index {output.bam}
		'''

rule shift_reads:
	input:
		bed_file=f'{data_dir}/{{cell_line}}/bed_combined/{{data_type}}.bed',
		phantompeakqualtools=f'{data_dir}/phantompeakqualtools.txt'
	output:
		f'{data_dir}/{{cell_line}}/bed_shifted/{{data_type}}.bed'
	run:
		# Input and DNase-seq are not shifted.
		if 'Input' in input.bed_file or 'Dnase' in input.bed_file:
			shell(f'cp {input.bed_file} {output}')
		else:
			# For MNase-seq data the shift is 149.
			if 'Nsome' in input.bed_file:
				shift = 149
			else:
				# Lookup the required shift produced by phantompeakqualtools.
				# These shifts are divided by 2 and rounded to the nearest integer.
				shifts = pd.read_csv(f'{data_dir}/phantompeakqualtools.txt', sep='\t', index_col=[0, 1], usecols=[0, 1, 3])
				shift = shifts.loc[(f'{wildcards.cell_line}', f'{wildcards.data_type}.bam')][0]
				shift = round(shift / 2)
			shell(f'shiftBed -i {input.bed_file} -g {genome_file} -p {shift} -m -{shift} > {output}')


# Step 5: Convert bed files to format accepted by RFECS
rule convert_bed_for_RFECS:
	input: f'{data_dir}/{{cell_line}}/bed_shifted/{{data_type}}.bed'
	output: f'{data_dir}/{{cell_line}}/bed_shifted_RFECS/{{data_type}}.bed'
	shell:
		r"""
		sed 's/ \+//g' {input} > {output}
		awk -F '\t' 'BEGIN {{OFS="\t"}} {{ if (($1) != "chrM")  print }}' {output} > {data_dir}/{wildcards.cell_line}/bed_shifted_RFECS/tmp_{wildcards.data_type}.bed
		mv {data_dir}/{wildcards.cell_line}/bed_shifted_RFECS/tmp_{wildcards.data_type}.bed {output}
		awk -F'\t' 'BEGIN{{OFS="\t"}}{{$5=""; gsub(FS"+",FS); print $0}}' {output} > {data_dir}/{wildcards.cell_line}/bed_shifted_RFECS/tmp_{wildcards.data_type}.bed
		mv {data_dir}/{wildcards.cell_line}/bed_shifted_RFECS/tmp_{wildcards.data_type}.bed {output}
		"""
