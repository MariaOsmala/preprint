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
threshold = 5

print(window)
print(bin_size)
print(N)
print(path)
print(cell_line)
print(normalizeBool)
print(NormCellLine)

source('code/find_random_with_signal.R')

# Loads profiles and regions
load(file = paste0(path, "/", cell_line, "/data_R/whole_genome_coverage2.RData"))

# Blacklist to remove problematic regions
DAC <- rtracklayer::import(paste0(path, "/blacklists/wgEncodeDacMapabilityConsensusExcludable.bed.gz"))
Duke <- rtracklayer::import(paste0(path, "/blacklists/wgEncodeDukeMapabilityRegionsExcludable.bed.gz"))
blacklist = union(DAC, Duke)
strand(blacklist) <- "+"

# Stay away from p300 peaks
p300 <- rtracklayer::import(paste0(path, '/K562/raw_data/wgEncodeAwgTfbsSydhK562P300IggrabUniPk.narrowPeak.gz'))

# Stay away from TSS annotations
TSS <- readRDS(paste0(path, "/GENCODE_TSS/", "GR_Gencode_protein_coding_TSS.RDS"))

print('Finding random with signal')
random_with_signal <- find_random_with_signal(profiles, regions, threshold,
                                              window, N, p300 = p300, TSS = TSS,
                                              blacklist = blacklist)
