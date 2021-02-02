# library(Rsamtools)
# library(snow)
# library(spp)
# library(accelerometry)
# library(Biostrings)
# library(bitops)
# library(BSgenome.Hsapiens.UCSC.hg19)
# library(circlize)
# library(doParallel)
# library(foreach)
# library(gdata)
# library(GenomicRanges)
# library(GetoptLong)
# library(ggplot2)
# library(grid)
# library(gridExtra)
# library(MASS)
library(optparse)
# library(pryr)
# library(RColorBrewer)
# library(reshape2)
# library(ROCR)
# library(rtracklayer)
# library(ShortRead)
# library(stringr)

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

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);


#if (is.null(opt)){
#  print_help(opt_parser)
#  stop("At least one argument must be supplied (input file).n", call.=FALSE)
#}


window=opt$window					#5000			
bin_size=opt$binSize      #100
N=opt$N                 #1000
path_to_dir=opt$pathToDir
cell_line=opt$cellLine
normalizeBool=opt$normalize
NormCellLine=opt$NormCellLine


print(window)				#5000			
print(bin_size)
print(N)                #1000
print(path_to_dir)
print(cell_line)
print(normalizeBool)
print(NormCellLine)

source("code/functions.R")
#setwd(path_to_dir)



#load the counts
path=path_to_dir
if(normalizeBool==TRUE){
	load( file=paste0(path,"/",NormCellLine,"/data_R/",N,"_enhancers_bin_",bin_size,"_window_",window,".RData")) #profiles, 	normalized_profiles, regions,
	countsOtherCellLine=profiles$counts #remove "V2" from the names
	names(countsOtherCellLine)=gsub("V2","",  names(countsOtherCellLine))

	load( file=paste0(path,"/",cell_line,"/data_R/",NormCellLine,"_normalized_",N,"_enhancers_bin_",bin_size,"_window_",window,".RData")) #profiles, normalized_profiles, regions,
	counts=profiles$counts #remove "V2" from the names
}else{
    load( file=paste0(path,"/",cell_line,"/data_R/",N,"_enhancers_bin_",bin_size,"_window_",window,".RData")) #profiles, normalized_profiles, regions,
    counts=profiles$counts #remove "V2" from the names
}
load(file=paste0(path,"/",cell_line,"/data_R/",N,"_randomRegions_with_signal_bin_",bin_size,"_window_",window,".RData")) #regions, accepted_GRanges,

################################3extract profiles#####################################################################3


strand(regions)="*"

figure_path=paste0(path_to_dir,"/figures/")

bam_folder=paste0(path,"/",cell_line,"/bam_shifted")
system.time(random_profiles<-extract_profiles_parallel(bam_folder=bam_folder, 
                                                                  regions=regions, directionality=FALSE, 
                                                                  directions=strand(regions), window=window, bin_size=bin_size))



bai_files=dir(pattern=".bai")
bam_files=dir(pattern=".bam")
bam_files=bam_files[which(bam_files %in% bai_files ==FALSE)]

setwd(path_to_dir)


mark_name<-NULL
for(i in 1:length(bam_files)){
  mark_name[i]=strsplit(bam_files[i],".bam")[[1]]
}

profiles<-list()
counts<-list()
profiles_temp<-list()
for(i in 1:length(random_profiles)){
  profiles_temp[[mark_name[i]]]=random_profiles[[i]][[1]]
  counts[[mark_name[i]]]=random_profiles[[i]][[2]]
}

profiles$profiles=profiles_temp
profiles$counts=counts

round_logic=FALSE #TRUE FALSE #11		
negatives_to_zero=FALSE #TRUE

normalized_profiles=normalize(profiles, round_logic, negatives_to_zero)

if(normalizeBool==TRUE){
  
  otherCellLine_normalized_profiles=normalize_Between_Cell_Lines(profiles=normalized_profiles, counts=counts, 
                                                        countsOtherCellLine=countsOtherCellLine, mods=names(normalized_profiles), 
                                                        round=FALSE, negatives_to_zero=FALSE)
  rm(counts, profiles_temp, random_profiles)
  
  
  #save counts and enhancer profiles
  
  print(N)
  
  save(profiles, normalized_profiles, otherCellLine_normalized_profiles, regions, accepted_GRanges, 
       file=paste0(path,"/",cell_line,"/data_R/",NormCellLine,"_normalized_",N,"_random_with_signal_bin_",bin_size,"_window_",window,".RData"))
  
  
  
}else{
  
  rm(counts, profiles_temp, random_profiles)
  #save counts and enhancer profiles
  
  print(N)
  
  save(profiles, normalized_profiles, regions, accepted_GRanges, file=paste0(path,"/",
                          cell_line,"/data_R/",N,"_random_with_signal_bin_",bin_size,"_window_",window,".RData"))
  
  
  
}


