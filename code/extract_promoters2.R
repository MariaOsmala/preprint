library(BSgenome.Hsapiens.UCSC.hg19)
library(optparse)

option_list = list(
  make_option(c("-w", "--window"), type="integer", default=5000, 
              help="window size [default=%default]", metavar="integer"),
  make_option(c("-b", "--binSize"), type="integer", default=100, help="bin size (resolution) [default= %default]", metavar="integer"),
  make_option(c("-N", "--N"), type="integer", default=1000000, 
              help="number of regions [default= %default]", metavar="integer"),
  make_option("--tssdist", type="integer", default=10000, 
              help="min distance between two TSS [default= %default]", metavar="integer"),
  make_option("--pathToDir", type="character", default="", 
              help="path to main folder [default= %default]", metavar="character"),
  make_option("--cellLine", type="character", default="", 
              help="cell line [default= %default]", metavar="character"),
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
# window=opt$window								
# bin_size=opt$binSize
# between_TSS_distance=opt$tssdist		
# N=opt$N 
# path_to_dir=opt$pathToDir
# cell_line=opt$cellLine
# DNase_peaks_file=opt$DNaseFile
# normalizeBool=opt$normalize
# NormCellLine=opt$NormCellLine

window=2000
bin_size=100
between_TSS_distance=2000
N=1000
path_to_dir='~/scratch_cs/csb/projects/enhancer_prediction/aaltorse/Data'
cell_line='K562'
DNase_peaks_file='~/scratch_cs/csb/projects/enhancer_prediction/aaltorse/Data/K562/raw_data/wgEncodeOpenChromDnaseK562PkV2.narrowPeak.gz'
NormCellLine=""

print(window)
print(bin_size)
print(between_TSS_distance)
print(N) 
print(path_to_dir)
print(cell_line)
print(DNase_peaks_file)
print(normalizeBool)
print(NormCellLine)


source('code/find_promoters.R')
source('code/create_profile.R')

directionality=TRUE # FALSE #


# Read promoter information
TSS_annotation <- readRDS(paste0(path, "/GENCODE_TSS/GR_Gencode_protein_coding_TSS.RDS"))
seqlengths(TSS_annotation) <- seqlengths(Hsapiens)[seqnames(seqinfo(TSS_annotation))]

# Used for finding promoter sites
DNase <- rtracklayer::import(DNase_peaks_file)

# Blacklist to remove problematic promoter sites
DAC <- rtracklayer::import(paste0(path, "/blacklists/wgEncodeDacMapabilityConsensusExcludable.bed.gz"))
Duke <- rtracklayer::import(paste0(path, "/blacklists/wgEncodeDukeMapabilityRegionsExcludable.bed.gz"))
blacklist = union(DAC, Duke)
strand(blacklist) <- "+"

promoters <- find_promoters(TSS_annotation = TSS_annotation,
                            DNase = DNase,
                            between_TSS_distance = between_TSS_distance, 
                            blacklist = blacklist,
                            window = window,
                            N = N)

# Create profiles for each promoter sites using all the histones available
profiles = list()

# First, collect a list of all the histone files
bam_files <- dir(paste0(path, '/', cell_line, '/bam_shifted'), pattern = "\\.bam$", full.name = TRUE)

# These are used to normalize the profiles
print('Reading control histone')
control_ind <- grep('Control', bam_files)
control <- BamFile(bam_files[control_ind])
profiles[['Control']] <- create_profile(promoters, histone = control, reference = NULL, ignore_strand = FALSE)

print('Reading input polymerase')
input_ind <- grep('Input', bam_files)
input <- BamFile(bam_files[input_ind])
profiles[['Input']] <- create_profile(promoters, histone = input, reference = NULL, ignore_strand = FALSE)

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
        histone <- BamFile(bam_file)
        profiles[[name]] <- create_profile(promoters, histone = histone,
                                           reference = reference, ignore_strand = FALSE)
    }
}

#normalize wrt to other cell line if applicaple
if(normalizeBool==TRUE){
    # Not implemented yet
}

# Save the profiles
dir.create(paste0(path, "/", cell_line, "/data_R"), recursive = TRUE, showWarnings = FALSE)
save(profiles, profiles_undirected, file = paste0(path, "/", cell_line,"/data_R/",N,"_promoters_bin_",bin_size,"_window_",window,"2.RData"))
