library(yaml)
library(argparser)
library(BSgenome.Hsapiens.UCSC.hg19)
library(GenomicRanges)
library(Rsamtools)

config <- read_yaml('workflow/config.yaml')

# parser <- arg_parser('Create profiles along the entire genome.')
# parser <- add_argument(parser, 'cell_line', type = 'character',
#                        help = paste0('The cell line to process. Either ', paste(config$cell_lines, collapse = ' or ')))
# cell_line <- parse_args(parser)$cell_line
cell_line <- 'K562'

source('code/profiles.R')
source('code/fname.R')

chr_names <- c('chr1','chr2','chr3',  'chr4',  'chr5','chr6',  'chr7',  'chr8',
               'chr9',  'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15',
               'chr16', 'chr17','chr18', 'chr19',   'chr20', 'chr21', 'chr22',
               'chrX') 
regions <- GRanges(chr_names, IRanges(start = 1, end = seqlengths(Hsapiens)[chr_names]))
regions <- slidingWindows(regions, width = opt$bin, step = opt$bin)
regions <- regions[width(regions) == opt$bin]  # Remove incomplete bins
regions <- unlist(regions)

# First, collect a list of all the BAM files
bam_files <- dir(opt$bam_folder, pattern = '\\.bam$', full.name = TRUE)

# These are used to normalize the profiles
cat('Computing whole genome coverage for control histone...\n')
control_ind <- grep('Control', bam_files)
control <- BamFile(bam_files[control_ind])
yieldSize(control) <- 1E6L
profiles_control <- create_profiles(regions, bam_file = control, reference = NULL, ignore_strand = TRUE)

cat('Computing whole genome coverage for input polymerase...\n')
input_ind <- grep('Input', bam_files)
input <- BamFile(bam_files[input_ind])
yieldSize(input) <- 1E6L
profiles_input <- create_profiles(regions, bam_file = input, reference = NULL, ignore_strand = TRUE)

# Create profiles for the rest of the BAM files
profiles = list()
for (bam_file in bam_files) {
    name <- tools::file_path_sans_ext(basename(bam_file))

    # Reference profiles have been created already, and we don't need Nsome
    # profiles.
    if (grepl('Input|Control|Nsome', name)) {
        next
    }

    # Determine reference profile
    if (grepl('Dnase', name)) {
        reference <- NULL
    } else if (grepl('Pol', name)) {
        reference <- profiles_input
    } else {
        reference <- profiles_control
    }

    # Create the profile
    cat(paste0('Computing whole genome coverage for: ', name, '\n'))
    bam_file <- BamFile(bam_file)
    yieldSize(bam_file) <- 1E6L
    profiles[[name]] <- create_profiles(regions, bam_file = bam_file,
                                        reference = reference, ignore_strand = TRUE)
}

# Make one big Profiles object
profiles <- do.call(cbind, profiles)

# Normalize with other cell line if applicable
if (!is.null(opt$normalize)) {
    reference <- readRDS(opt$normalize)
    profiles <- normalize_profiles(profiles, reference)
}

saveRDS(profiles, file = opt$output)
cat(paste0('Saved whole genome coverage data to ', opt$output, '\n'))
