library(BSgenome.Hsapiens.UCSC.hg19)
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

source('code/create_profile.R')

#############################Sample random regions###################################################
#compute the probability for each chromosome, include only chr1-21 and chrX

allowed_chroms=c("chr1",  "chr2",  "chr3",  "chr4",  "chr5",  "chr6",  "chr7",  "chr8",  "chr9", 
"chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18",
"chr19", "chr20", "chr21", "chr22", "chrX") 

whole_genome <- GRanges(allowed_chroms, IRanges(start=1, end=seqlengths(Hsapiens)[allowed_chroms]))

# Add the blacklist to the mask of regions not to use
DAC <- rtracklayer::import(paste0(path, "/blacklists/wgEncodeDacMapabilityConsensusExcludable.bed.gz"))
Duke <- rtracklayer::import(paste0(path, "/blacklists/wgEncodeDukeMapabilityRegionsExcludable.bed.gz"))
blacklist <- union(DAC, Duke)

# Add p300 BSs to the mask of regions not to use
p300 <- rtracklayer::import(p300_peaks_file)
start(p300) <- start(p300) + p300$peak
p300 <- resize(p300, 1000, fix="center")
blacklist <- union(blacklist, p300)

# Add protein coding TSS to the mask of regions not to use
TSS <- readRDS(paste0(path, "/GENCODE_TSS/", "GR_Gencode_protein_coding_TSS.RDS"))
blacklist <- union(blacklist, TSS)

blacklist <- keepSeqlevels(blacklist, allowed_chroms, pruning.mode = "coarse")
random_regions <- regioneR::createRandomRegions(N, length.mean = window, length.sd = 0, genome = whole_genome, mask = blacklist)

# Create profiles for each region using all the histones available
profiles = list()
profiles_undirected = list()

# First, collect a list of all the histone files
bam_files <- dir(paste0(path, '/', cell_line, '/bam_shifted'), pattern = "\\.bam$", full.name = TRUE)

# These are used to normalize the profiles
print('Reading control histone')
control_ind <- grep('Control', bam_files)
control <- rtracklayer::import(bam_files[control_ind])
profiles[['Control']] <- create_profile(random_regions, histone = control, reference = NULL, ignore_strand = FALSE)

print('Reading input polymerase')
input_ind <- grep('Input', bam_files)
input <- rtracklayer::import(bam_files[input_ind])
profiles[['Input']] <- create_profile(random_regions, histone = input, reference = NULL, ignore_strand = FALSE)

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
        histone <- rtracklayer::import(bam_file)
        profiles[[name]] <- create_profile(random_regions, histone = histone,
                                           reference = reference, ignore_strand = FALSE)
    }
}

#normalize wrt to other cell line if applicaple
if(normalizeBool==TRUE){
    # Not implemented yet
}

# Save the profiles
dir.create(paste0(path, "/", cell_line, "/data_R"), recursive = TRUE, showWarnings = FALSE)
save(profiles, profiles_undirected, file = paste0(path, "/", cell_line,"/data_R/",N,"_pure_random_1000_bin_",bin_size,"_window_",window,"2.RData"))
