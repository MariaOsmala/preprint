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
# 
# 
# #SAVE counts
# if(N==1000000){
#   N=length(promoters$promoters)
# }
# print("number of promoters")
# print(N)
# 
# promoters=promoters[1:N]
# #promoters$strand=promoters$strand[1:N]
# 
# strand(promoters)="*" #THIS IS IMPORTANT
# 
# promoter_profiles_directed<-extract_profiles_parallel(bam_folder, regions=promoters, directionality=TRUE, 
#                                                       directions=promoters$direction, window, bin_size)
# promoter_profiles_undirected<-extract_profiles_parallel(bam_folder, regions=promoters, directionality=FALSE, 
#                                                         directions=promoters$direction, window, bin_size)
# 
# setwd(bam_folder)
# bai_files=dir(pattern=".bai")
# bam_files=dir(pattern=".bam")
# bam_files=bam_files[which(bam_files %in% bai_files ==FALSE)]
# 
# setwd(path_to_dir)
# 
# 
# mark_name<-NULL
# 
# for(i in 1:length(bam_files)){
#   mark_name[i]=strsplit(bam_files[i],".bam")[[1]]
# }
# 
# profiles_directed<-list()
# counts<-list()
# profiles_temp<-list()
# for(i in 1:length(promoter_profiles_directed)){
#   profiles_temp[[mark_name[i]]]=promoter_profiles_directed[[i]][[1]]
#   counts[[mark_name[i]]]=promoter_profiles_directed[[i]][[2]]
# }
# 
# profiles_directed$profiles=profiles_temp
# profiles_directed$counts=counts
# 
# profiles_undirected<-list()
# counts<-list()
# profiles_temp<-list()
# for(i in 1:length(promoter_profiles_undirected)){
#   profiles_temp[[mark_name[i]]]=promoter_profiles_undirected[[i]][[1]]
#   counts[[mark_name[i]]]=promoter_profiles_undirected[[i]][[2]]
# }
# 
# profiles_undirected$profiles=profiles_temp
# profiles_undirected$counts=counts
# 
# normalized_profiles_directed=normalize(profiles_directed, round_logic, negatives_to_zero)
# normalized_profiles_undirected=normalize(profiles_undirected, round_logic, negatives_to_zero)
# 
# if(normalizeBool==TRUE){
#   otherCellLine_normalized_profiles_directed=normalize_Between_Cell_Lines(profiles=normalized_profiles_directed, counts=counts, 
#                                                                  countsOtherCellLine=countsOtherCellLine, mods=names(normalized_profiles_directed), 
#                                                                  round=FALSE, negatives_to_zero=FALSE)
#   otherCellLine_normalized_profiles_undirected=normalize_Between_Cell_Lines(profiles=normalized_profiles_undirected, counts=counts, 
#                                                                    countsOtherCellLine=countsOtherCellLine, mods=names(normalized_profiles_undirected), 
#                                                                    round=FALSE, negatives_to_zero=FALSE)
#   
#   rm(counts, profiles_temp, promoter_profiles_directed, promoter_profiles_undirected)
#   
#   
#   #save counts and enhancer profiles
#   regions=promoters
#   rm(promoters)
#   
#   save(profiles_directed,profiles_undirected, normalized_profiles_directed,normalized_profiles_undirected,
#        otherCellLine_normalized_profiles_directed,otherCellLine_normalized_profiles_undirected, 
#        regions, file=paste0(path, "/", cell_line, "/data_R/", NormCellLine, "_normalized_", N, "_promoters_bin_", bin_size, "_window_", window, ".RData"))
#   
#   
#   
# }else{
#   rm(counts, profiles_temp, promoter_profiles_directed, promoter_profiles_undirected)
#   
#   
#   #save counts and enhancer profiles
#   regions=promoters
#   rm(promoters)
#   
#   save(profiles_directed,profiles_undirected, normalized_profiles_directed,normalized_profiles_undirected, 
#        regions, file=paste0(path, "/", cell_line, "/data_R/", N, "_promoters_bin_", bin_size, "_window_", window, ".RData"))
#   
#   
# }
# 











