library(BSgenome.Hsapiens.UCSC.hg19)
library(GenomicRanges)
library(Rsamtools)
library(optparse)

option_list = list(
  make_option(c("-w", "--window"), type="integer", default=5000, 
              help="window size [default=%default]", metavar="integer"),
  make_option(c("-b", "--binSize"), type="integer", default=100, help="bin size (resolution) [default= %default]", metavar="integer"),
  make_option(c("-N", "--N"), type="integer", default=1000000, 
              help="number of regions [default= %default]", metavar="integer"),
  make_option("--pathToDir", type="character", default="", 
              help="path to main folder [default= %default]", metavar="character"),
  make_option("--cellLine", type="character", default="", 
              help="cell line [default= %default]", metavar="character"),
  make_option("--normalize", type="logical", default=FALSE, 
              help="do we normalize wrt data from other cell line [default= %default]", metavar="logical"),
  make_option("--NormCellLine", type="character", default="", 
              help="name of the cell line normalized wrt [default= %default]", metavar="character")
); 

# opt_parser = OptionParser(option_list=option_list);
# opt = parse_args(opt_parser);
# 
# 
# window=opt$window					#5000			
# bin_size=opt$binSize      #100
# N=opt$N                 #1000
# path_to_dir=opt$pathToDir
# cell_line=opt$cellLine
# normalizeBool=opt$normalize
# NormCellLine=opt$NormCellLine

window=2000
bin_size=100
N=1000
path='~/scratch_cs/csb/projects/enhancer_prediction/aaltorse/Data'
cell_line='K562'
normalizeBool=FALSE
NormCellLine=NULL

print(window)
print(bin_size)
print(N)
print(path)
print(cell_line)
print(normalizeBool)
print(NormCellLine)


chr_names <- c("chr1","chr2","chr3",  "chr4",  "chr5","chr6",  "chr7",  "chr8",
               "chr9",  "chr10", "chr11", "chr12", "chr13", "chr14", "chr15",
               "chr16", "chr17","chr18", "chr19",   "chr20", "chr21", "chr22",
               "chrX") 
regions <- GRanges(chr_names, IRanges(start=1, end=seqlengths(Hsapiens)[chr_names]))
regions <- slidingWindows(regions, width = bin_size, step = bin_size)
regions <- regions[width(regions) == bin_size]  # Remove incomplete bins
regions <- unlist(regions)


source('code/create_profiles.R')

# First, collect a list of all the BAM files
bam_files <- dir(paste0(path, '/', cell_line, '/bam_shifted'), pattern = "\\.bam$", full.name = TRUE)
# bam_files <- dir(paste0('~/scratch_cs/csb/projects/enhancer_prediction/experiments/RProjects/preprint/Data', '/', cell_line, '/bam_shifted'), pattern = "\\.bam$", full.name = TRUE)

# These are used to normalize the profiles
print('Computing whole genome coverage for control histone')
control_ind <- grep('Control', bam_files)
control <- BamFile(bam_files[control_ind])
yieldSize(control) <- 1E6L
profiles_control <- create_profiles(regions, bam_file = control, reference = NULL, ignore_strand = TRUE)

print('Computing whole genome coverage for input polymerase')
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
    print(paste0("Computing whole genome coverage for: ", name))
    bam_file <- BamFile(bam_file)
    yieldSize(bam_file) <- 1E6L
    profiles[[name]] <- create_profiles(regions, bam_file = bam_file,
                                        reference = reference, ignore_strand = TRUE)
}

#normalize wrt to other cell line if applicable
if(normalizeBool==TRUE){
    # Not implemented yet
}

save(profiles, file = paste0(path, "/", cell_line, "/data_R/whole_genome_coverage.RData"))
