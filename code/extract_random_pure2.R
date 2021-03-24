library(Rsamtools)
library(optparse)

# option_list = list(
#   make_option(c("-w", "--window"), type="integer", default=5000, 
#               help="window size [default=%default]", metavar="integer"),
#   make_option(c("-b", "--binSize"), type="integer", default=100, help="bin size (resolution) [default= %default]", metavar="integer"),
#   make_option(c("-N", "--N"), type="integer", default=1000000, 
#               help="number of regions [default= %default]", metavar="integer"), 
#   make_option("--pathToDir", type="character", default="", 
#               help="path to main folder [default= %default]", metavar="character"),
#   make_option("--p300File", type="character", default="", 
#               help="path to p300 peak file [default= %default]", metavar="character"),
#   make_option("--cellLine", type="character", default="", 
#               help="cell line [default= %default]", metavar="character"),
#   make_option("--normalize", type="logical", default=FALSE, 
#               help="do we normalize wrt data from other cell line [default= %default]", metavar="logical"),
#   make_option("--NormCellLine", type="character", default="", 
#               help="name of the cell line normalized wrt [default= %default]", metavar="character")
# ); 
# 
# opt_parser = OptionParser(option_list=option_list);
# opt = parse_args(opt_parser);

#if (is.null(opt)){
#  print_help(opt_parser)
#  stop("At least one argument must be supplied (input file).n", call.=FALSE)
#}


#window=opt$window					#5000			
#bin_size=opt$binSize      #100
#N=opt$N
#path_to_dir=opt$pathToDir
#cell_line=opt$cellLine
#p300_peaks_file=opt$p300File
#normalizeBool=opt$normalize
#NormCellLine=opt$NormCellLine

window=2000
bin_size=100
N=1000
path_to_dir='~/scratch_cs/csb/projects/enhancer_prediction/aaltorse/Data'
cell_line='K562'
p300_peaks_file='~/scratch_cs/csb/projects/enhancer_prediction/aaltorse/Data/K562/raw_data/wgEncodeAwgTfbsSydhK562P300IggrabUniPk.narrowPeak.gz'
normalizeBool=FALSE
NormCellLine=NULL
max_dist_to_promoter=2000

source('code/find_random.R')
source('code/create_profile.R')

# Add the blacklist to the mask of regions not to use
DAC <- rtracklayer::import(paste0(path, "/blacklists/wgEncodeDacMapabilityConsensusExcludable.bed.gz"))
Duke <- rtracklayer::import(paste0(path, "/blacklists/wgEncodeDukeMapabilityRegionsExcludable.bed.gz"))
blacklist <- union(DAC, Duke)

# Stay away from p300 peaks
p300 <- rtracklayer::import(paste0(path, '/K562/raw_data/wgEncodeAwgTfbsSydhK562P300IggrabUniPk.narrowPeak.gz'))

# Stay away from TSS annotations
TSS <- readRDS(paste0(path, "/GENCODE_TSS/", "GR_Gencode_protein_coding_TSS.RDS"))

print('Finding random pure')
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
save(profiles, file = paste0(path, "/", cell_line,"/data_R/",N,"_pure_random_1000_bin_",bin_size,"_window_",window,"2.RData"))
