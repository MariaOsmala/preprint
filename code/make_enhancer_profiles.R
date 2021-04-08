library(GenomicRanges)
library(Rsamtools)

path <- '~/scratch_cs/csb/projects/enhancer_prediction/aaltorse/Data'
# cell_line <- 'K562'
cell_line <- 'Gm12878'
window <- 2000
N <- 1000
max_dist_to_promoter <- 2000
between_TSS_distance <- 2000
bin_size <- 100
threshold <- 5


source('code/find_enhancers.R')
source('code/find_promoters.R')
source('code/find_random.R')
source('code/TSS_protein_coding.R')
source('code/profiles.R')

# These are the chromosomes of interest
chroms_of_interest = c('chr1',  'chr2',  'chr3',  'chr4',  'chr5',  'chr6',
                       'chr7',  'chr8',  'chr9', 'chr10', 'chr11', 'chr12',
                       'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18',
                       'chr19', 'chr20', 'chr21', 'chr22', 'chrX') 

# Read TSS annotations to extract protein codings
TSS_annotation <- TSS_protein_coding(paste0(path, '/GENCODE_TSS/gencode.v27lift37.annotation.gtf.gz'))

cat('Reading p300 and DNase peaks...')
p300 <- rtracklayer::import(paste0(path, '/', cell_line, '/raw_data/wgEncodeAwgTfbsSydhK562P300IggrabUniPk.narrowPeak.gz'))
DNase <- rtracklayer::import( paste0(path, '/', cell_line, '/raw_data/wgEncodeOpenChromDnaseK562PkV2.narrowPeak.gz'))
cat(' done.\n')

# Blacklist to remove problematic sites
cat('Creating blacklist...')
DAC <- rtracklayer::import(paste0(path, '/blacklists/wgEncodeDacMapabilityConsensusExcludable.bed.gz'))
Duke <- rtracklayer::import(paste0(path, '/blacklists/wgEncodeDukeMapabilityRegionsExcludable.bed.gz'))
blacklist = union(DAC, Duke)
cat(' done.\n')

# Find interesting sites
enhancers <- find_enhancers(p300, DNase, window = window, N = N,
                            TSS = TSS_annotation, max_dist_to_promoter = max_dist_to_promoter,
                            blacklist = blacklist)

promoters <- find_promoters(TSS_annotation, DNase, 
                            between_TSS_distance = between_TSS_distance, 
                            window = window, N = N,
                            blacklist = blacklist)

random_regions <- find_random(window, N, p300 = p300, TSS = TSS_annotation,
                              chroms_of_interest = chroms_of_interest, blacklist = blacklist)

# Load whole genome coverage
load(file = paste0(path, '/', cell_line, '/data_R/whole_genome_coverage2.RData'))
coverage <- do.call(cbind, profiles)  # Concatenate to one big matrix

# Blacklist all regions which don't have enough coverage
coverage[coverage < 0] <- 0
blacklist <- union(blacklist, regions[rowSums(coverage) < threshold])

random_regions_with_signal <- find_random(window, N, p300 = p300, TSS = TSS_annotation,
                                          chroms_of_interest = chroms_of_interest, blacklist = blacklist)

# All sites
sites <- c(enhancers, promoters, random_regions, random_regions_with_signal)

# Create profiles for each site using all the histones available
profiles = list()

# First, collect a list of all the histone files
bam_files <- dir(paste0(path, '/', cell_line, '/bam_shifted'), pattern = '\\.bam$', full.name = TRUE)

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
for (bam_file in bam_files) {
    name <- tools::file_path_sans_ext(basename(bam_file))

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
        profiles[[name]] <- create_profiles(sites, bam_file = BamFile(bam_file),
                                            reference = reference, ignore_strand = TRUE)
        cat(' done.\n')
    }
}

#normalize wrt to other cell line if applicaple
# if(normalizeBool==TRUE){
#     # Not implemented yet
# }

# Make one big Profiles object
profiles <- do.call(cbind, profiles)

# Save the profiles
dir.create(paste0(path, '/', cell_line, '/data_R'), recursive = TRUE, showWarnings = FALSE)
profiles_file <- paste0(path, '/', cell_line, '/data_R/profiles.RData')
save(profiles, file = profiles_file)
cat(paste0('Data saved to ', profiles_file, '\n'))
