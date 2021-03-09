library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg19)
library(optparse)

option_list = list(
  make_option(c("-w", "--window"), type="integer", default=5000, 
              help="window size [default=%default]", metavar="integer"),
  make_option(c("-b", "--binSize"), type="integer", default=100, help="bin size (resolution) [default= %default]", metavar="integer"),
  make_option(c("-N", "--N"), type="integer", default=1000000, 
              help="number of regions [default= %default]", metavar="integer"),
  make_option(c("-p", "--distToPromoter"), type="integer", default=2000, 
              help="min distance to any TSS as bp [default= %default]", metavar="integer"),
  make_option("--pathToDir", type="character", default="", 
              help="path to main folder [default= %default]", metavar="character"),
  make_option("--cellLine", type="character", default="", 
              help="cell line [default= %default]", metavar="character"),
  make_option("--p300File", type="character", default="", 
              help="path to p300 peak file [default= %default]", metavar="character"),
  make_option("--DNaseFile", type="character", default="", 
              help="path to DNase peak file [default= %default]", metavar="character"),
  make_option("--normalize", type="logical", default=FALSE, 
              help="do we normalize wrt data from other cell line [default= %default]", metavar="logical"),
  make_option("--NormCellLine", type="character", default="", 
              help="name of the cell line normalized wrt [default= %default]", metavar="character")
  
); 

# opt_parser = OptionParser(option_list=option_list);
# opt = parse_args(opt_parser);
# 
# if (is.null(opt)){
#   print_help(opt_parser)
#   stop("At least one argument must be supplied (input file).n", call.=FALSE)
# }

# window=opt$window								
# distance_to_promoters=opt$distToPromoter										
# bin_size=opt$binSize
# N=opt$N
# path=opt$pathToDir
# cell_line=opt$cellLine
# p300_peaks_file=opt$p300File
# DNase_peaks_file=opt$DNaseFile
# 
# normalizeBool=opt$normalize
# 
# NormCellLine=opt$NormCellLine

window=2000
distance_to_promoters=2000
bin_size=100
N=1000
path='~/scratch_cs/csb/projects/enhancer_prediction/aaltorse/Data'
cell_line='K562'
p300_peaks_file='~/scratch_cs/csb/projects/enhancer_prediction/aaltorse/Data/K562/raw_data/wgEncodeAwgTfbsSydhK562P300IggrabUniPk.narrowPeak.gz'
DNase_peaks_file='~/scratch_cs/csb/projects/enhancer_prediction/aaltorse/Data/K562/raw_data/wgEncodeOpenChromDnaseK562PkV2.narrowPeak.gz'
normalizeBool=FALSE
NormCellLine=""

print(window)
print(distance_to_promoters)
print(bin_size)
print(N)
print(path)
print(cell_line)
print(p300_peaks_file)
print(DNase_peaks_file)
print(normalizeBool)
print(NormCellLine)

source('code/find_enhancers.R')
source('code/create_profile.R')

# Used for finding enhancer sites
p300 <- rtracklayer::import(p300_peaks_file)
DNase <- rtracklayer::import(DNase_peaks_file)

# Promoters used for eliminating enhancer sites that are too close
promoters <- readRDS(paste0(path, "/GENCODE_TSS/GR_Gencode_protein_coding_TSS_positive.RDS"))
seqlengths(promoters) <- seqlengths(Hsapiens)[seqnames(seqinfo(promoters))]

# Blacklist to remove problematic enhancer sites
DAC <- rtracklayer::import(paste0(path, "/blacklists/wgEncodeDacMapabilityConsensusExcludable.bed.gz"))
Duke <- rtracklayer::import(paste0(path, "/blacklists/wgEncodeDukeMapabilityRegionsExcludable.bed.gz"))
blacklist = union(DAC, Duke)
strand(blacklist) <- "+"

# Find the enhancer sites
enhancers  <- find_enhancers(p300, DNase, window = window, N = N,
                             promoters = promoters, max_dist_to_promoter = distance_to_promoters,
                             blacklist = blacklist)

# Create profiles for each enhancer sites using all the histones available
profiles = list()

# First, collect a list of all the histone files
bam_files <- dir(paste0(path, '/', cell_line, '/bam_shifted'), pattern = "\\.bam$", full.name = TRUE)

# These are used to normalize the profiles
print('Reading control histone')
control_ind <- grep('Control', bam_files)
profiles[['Control']] <- create_profile(enhancers, histone = rtracklayer::import(bam_files[control_ind]),
                                        reference = NULL)
print('Reading input polymerase')
input_ind <- grep('Input', bam_files)
profiles[['Input']] <- create_profile(enhancers, histone = rtracklayer::import(bam_files[input_ind]),
                                      reference = NULL)

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
        profiles[[name]] <- create_profile(enhancers, histone = rtracklayer::import(bam_file),
                                           reference = reference)
    }
}

#normalize wrt to other cell line if applicaple
if(normalizeBool==TRUE){
    # Not implemented yet
}

# Save the profiles
# dir.create(paste0(path, "/", cell_line, "/data_R"), recursive = TRUE, showWarnings = FALSE)
# save(profiles, file = paste0(path, "/", cell_line,"/data_R/",N,"_enhancers_bin_",bin_size,"_window_",window,".RData"))
