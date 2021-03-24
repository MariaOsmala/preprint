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
threshold = 5

print(window)
print(bin_size)
print(N)
print(path)
print(cell_line)
print(normalizeBool)
print(NormCellLine)

source('code/find_random.R')
source('code/create_profile.R')

# Loads profiles and regions
load(file = paste0(path, "/", cell_line, "/data_R/whole_genome_coverage2.RData"))

# Blacklist to remove problematic regions
DAC <- rtracklayer::import(paste0(path, "/blacklists/wgEncodeDacMapabilityConsensusExcludable.bed.gz"))
Duke <- rtracklayer::import(paste0(path, "/blacklists/wgEncodeDukeMapabilityRegionsExcludable.bed.gz"))
blacklist = union(DAC, Duke)

# Select the $data parts of the coverage list
coverage <- unlist(profiles, recursive = FALSE)[c(TRUE, FALSE)]

# Concatenate to one big matrix
coverage <- t(do.call(rbind, coverage))

# Don't count negative coverage values 
coverage[coverage < 0] <- 0

# Blacklist all regions which don't have enough coverage
blacklist <- union(blacklist, regions[rowSums(coverage) < threshold])

# Stay away from p300 peaks
p300 <- rtracklayer::import(paste0(path, '/K562/raw_data/wgEncodeAwgTfbsSydhK562P300IggrabUniPk.narrowPeak.gz'))

# Stay away from TSS annotations
TSS <- readRDS(paste0(path, "/GENCODE_TSS/", "GR_Gencode_protein_coding_TSS.RDS"))

print('Finding random with signal')
random_regions <- find_random(window, N, p300 = p300, TSS = TSS,
                              blacklist = blacklist)

# Create profiles for each region using all the histones available
profiles = list()

# First, collect a list of all the histone files
bam_files <- dir(paste0(path, '/', cell_line, '/bam_shifted'), pattern = "\\.bam$", full.name = TRUE)

# These are used to normalize the profiles
print('Reading control histone')
control_ind <- grep('Control', bam_files)
profiles[['Control']] <- create_profile(random_regions, bam_file = BamFile(bam_files[control_ind]),
                                        reference = NULL, ignore_strand = TRUE)
print('Reading input polymerase')
input_ind <- grep('Input', bam_files)
profiles[['Input']] <- create_profile(random_regions, bam_file = BamFile(bam_files[input_ind]),
                                      reference = NULL, ignore_strand = TRUE)

# Create profiles for the rest of the histones
for (bam_file in bam_files) {
    name <- tools::file_path_sans_ext(basename(bam_file))

    # Determine reference profile
    if (length(grep('Dnase|Nsome', name)) > 0) {
        reference <- NULL
    } else if (length(grep('Pol', name)) > 0) {
        reference <- profiles[['Input']]
    } else {
        reference <- profiles[['Control']]
    }

    # Create the profile (reference profiles have been created already)
    if (length(grep('Input|Control', name)) == 0) {
        print(paste0("Processing: ", name))
        profiles[[name]] <- create_profile(random_regions, bam_file = BamFile(bam_file),
                                           reference = reference, ignore_strand = TRUE)
    }
}

#normalize wrt to other cell line if applicaple
if(normalizeBool==TRUE){
    # Not implemented yet
}

# Save the profiles
dir.create(paste0(path, "/", cell_line, "/data_R"), recursive = TRUE, showWarnings = FALSE)
save(profiles, file = paste0(path, "/", cell_line,"/data_R/",N,"_random_with_signal_1000_bin_",bin_size,"_window_",window,"2.RData"))
