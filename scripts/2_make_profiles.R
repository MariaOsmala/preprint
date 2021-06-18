library(GenomicRanges)
library(Rsamtools)
library(argparser)
library(yaml)
library(preprint)

source('R/find_enhancers.R')
config <- read_yaml('workflow/config.yaml')

# parser <- arg_parser('Create the profiles that will serve as training data.')
# parser <- add_argument(parser, 'cell_line', type = 'character',
#                        help = paste0('The cell line to process. Either ', paste(config$cell_lines, collapse = ' or ')))
# cell_line <- parse_args(parser)$cell_line
cell_line <- 'K562'

# These are the chromosomes of interest
chroms_of_interest = c('chr1',  'chr2',  'chr3',  'chr4',  'chr5',  'chr6',
                       'chr7',  'chr8',  'chr9', 'chr10', 'chr11', 'chr12',
                       'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18',
                       'chr19', 'chr20', 'chr21', 'chr22', 'chrX') 

# Read TSS annotations to extract protein codings
TSS_annotations <- TSS_protein_coding(paste0(config$data_dir, '/GENCODE_TSS/gencode.v27lift37.annotation.gtf.gz'))

cat('Reading p300 and DNase peaks...')
p300 <- rtracklayer::import(paste0(config$data_dir, '/', cell_line, '/raw_data/wgEncodeAwgTfbsSydhK562P300IggrabUniPk.narrowPeak.gz'))
DNase <- rtracklayer::import(paste0(config$data_dir, '/', cell_line, '/raw_data/wgEncodeOpenChromDnaseK562PkV2.narrowPeak.gz'))
cat(' done.\n')

# Blacklist to remove problematic sites
cat('Creating blacklist...')
blacklist = GRanges()
for (blacklist_type in c('DacMapabilityConsensus', 'DukeMapabilityRegions')) {
    blacklist_file <- paste0(config$data_dir, '/blacklists/wgEncode', blacklist_type, 'Excludable.bed.gz')
    blacklist = union(blacklist, rtracklayer::import(blacklist_file))
}
cat(' done.\n')

# Find interesting sites
enhancers <- find_enhancers(p300, DNase, window = config$profiles$window_size, N = config$profiles$num_enhancers,
                            TSS_annotations = TSS_annotations, min_dist_to_promoter = config$profiles$min_dist_to_promoter,
                            blacklist = blacklist)

promoters <- find_promoters(TSS_annotations, DNase, 
                            between_promoter_distance = config$profiles$min_dist_between_promoters,
                            window = config$profiles$window_size, N = config$profiles$num_promoters,
                            blacklist = blacklist)

random_regions <- find_random(config$profiles$window_size, N = config$profiles$num_random_pure,
                              p300 = p300, TSS = TSS_annotations,
                              chroms_of_interest = chroms_of_interest,
                              blacklist = blacklist)

# Load whole genome coverage
coverage <- readRDS(paste0(config$data_dir, '/', cell_line, '/data_R/whole_genome_coverage.rds'))

# Blacklist all regions which don't have enough coverage
not_enough_coverage <- rowSums(coverage) < config$profiles$signal_threshold
blacklist <- union(blacklist, attr(coverage, 'ranges')[not_enough_coverage])

# The whole genome coverage uses a lot of memory and is no longer needed.
# Clearing it now makes space for the profile data later, which is also memory hungry.
rm(coverage)

random_regions_with_signal <- find_random(config$profiles$window_size,
                                          config$profiles$num_random_with_signal,
                                          p300 = p300, TSS = TSS_annotations,
                                          chroms_of_interest = chroms_of_interest,
                                          blacklist = blacklist)

# All sites
sites <- c(enhancers, promoters, random_regions, random_regions_with_signal)

# Create profiles for each site using all the histones available
profiles = list()

# First, collect a list of all the histone files
bam_files <- dir(paste0(config$data_dir, '/', cell_line, '/bam_shifted'), pattern = '\\.bam$', full.name = TRUE)

# These are used to normalize the profiles
cat('Reading control histone...')
control_ind <- grep('Control', bam_files)
profiles_control <- create_profiles(sites, bam_file = BamFile(bam_files[control_ind]),
                                    reference = NULL, ignore_strand = TRUE)
cat(' done.\n')

cat('Reading input polymerase...')
input_ind <- grep('Input', bam_files)
profiles_input <- create_profiles(sites, bam_file = BamFile(bam_files[input_ind]),
                                  reference = NULL, ignore_strand = TRUE)
cat(' done.\n')

# Create profiles for the rest of the histones
for (filename in bam_files) {
    name <- tools::file_path_sans_ext(basename(filename))

    # Determine reference profile
    if (length(grep('Dnase|Nsome', name)) > 0) {
        reference <- NULL
    } else if (length(grep('Pol', name)) > 0) {
        reference <- profiles_input
    } else {
        reference <- profiles_control
    }

    # Create the profile (reference profiles have been created already)
    if (length(grep('Input|Control', name)) == 0) {
        cat(paste0('Processing ', name, '...'))
        profiles[[name]] <- create_profiles(sites, bam_file = BamFile(filename),
                                            reference = reference, ignore_strand = TRUE)
        cat(' done.\n')
    }
}

# Make one big Profiles object
profiles <- do.call(cbind, profiles)

# Save the profiles
fname <- paste0(config$data_dir, '/', cell_line, '/data_R/profiles.rds')
dir.create(dirname(fname), recursive = TRUE, showWarnings = FALSE)
saveRDS(profiles, file = fname)
cat(paste0('Data saved to ', fname, '\n'))
