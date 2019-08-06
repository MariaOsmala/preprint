library(Rsamtools)
library(snow)
library(spp)
library(accelerometry)
library(Biostrings)
library(bitops)
library(BSgenome.Hsapiens.UCSC.hg19)
library(circlize)
library(doParallel)
library(foreach)
library(gdata)
library(GenomicRanges)
library(GetoptLong)
library(ggplot2)
library(grid)
library(gridExtra)
library(MASS)
library(optparse)
library(pryr)
library(RColorBrewer)
library(reshape2)
library(ROCR)
library(rtracklayer)
library(ShortRead)
library(stringr)

option_list = list(
  make_option(c("-w", "--window"), type="integer", default=5000, 
              help="window size [default=%default]", metavar="integer"),
  make_option(c("-b", "--binSize"), type="integer", default=100, help="bin size (resolution) [default= %default]", metavar="integer"),
  make_option(c("-N", "--N"), type="integer", default=1000000, 
              help="number of regions [default= %default]", metavar="integer"),
  make_option(c("-p", "--distToPromoter"), type="integer", default=2000, 
              help="min distance to any TSS as bp [default= %default]", metavar="integer"),
  make_option(c("-pathToDir", "--pathToDir"), type="character", default="", 
              help="path to main folder [default= %default]", metavar="character"),
  make_option(c("-cellLine", "--cellLine"), type="character", default="", 
              help="cell line [default= %default]", metavar="character"),
  make_option(c("-p300File", "--p300File"), type="character", default="", 
              help="path to p300 peak file [default= %default]", metavar="character"),
  make_option(c("-DNaseFile", "--DNaseFile"), type="character", default="", 
              help="path to DNase peak file [default= %default]", metavar="character"),
  make_option(c("-normalize", "--normalize"), type="logical", default=FALSE, 
              help="do we normalize wrt data from other cell line [default= %default]", metavar="logical"),
  make_option(c("-NormCellLine", "--NormCellLine"), type="character", default="", 
              help="name of the cell line normalized wrt [default= %default]", metavar="character")
  
  
  
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
}

window=opt$window								
distance_to_promoters=opt$distToPromoter										
bin_size=opt$binSize
N=opt$N
path_to_dir=opt$pathToDir
cell_line=opt$cellLine
p300_peaks_file=opt$p300File
DNase_peaks_file=opt$DNaseFile

normalizeBool=opt$normalize

NormCellLine=opt$NormCellLine

setwd(path_to_dir)
source("code/functions.R")
source("code/extract_profiles_parallel.R")


path=paste(path_to_dir, "/Data/",sep="")
figure_path=paste(path_to_dir,"/figures/",sep="")

#TSS_annotation:
#GFF/GTF File format, start and end are both one-based

GR_Gencode_protein_coding_TSS_positive=readRDS(paste(path,"GENCODE_TSS/","GR_Gencode_protein_coding_TSS_positive.RDS",sep=""))

#these need the chromosome length information

human.chromlens = seqlengths(Hsapiens)

seqlengths(GR_Gencode_protein_coding_TSS_positive)<-human.chromlens[seqnames(seqinfo(GR_Gencode_protein_coding_TSS_positive))]

bam_folder=paste(path,cell_line,"/bam_shifted" ,sep="")

round_logic=FALSE #TRUE FALSE #11		
negatives_to_zero=FALSE #TRUE

#are we normalizing the data wrt some other cell line data
if(normalizeBool==TRUE){
  load( file=paste(path,NormCellLine,"/data_R/",N,"_enhancers_bin_",bin_size,"_window_",window,".RData",sep="")) #profiles, normalized_profiles, regions,
  countsOtherCellLine=profiles$counts #remove "V2" from the names
  names(countsOtherCellLine)=gsub("V2","",  names(countsOtherCellLine))
}

#-------------extract the training data enhancers from a file --------------------


ENCODE_blacklist=ENCODE_blaclist_regions(path_to_dir)
TSS_annotation=GR_Gencode_protein_coding_TSS_positive
enhancers <- create_enhancer_list(p300=p300_peaks_file, DNase_peaks_file=DNase_peaks_file,  #15488
                                  TSS_annotation=TSS_annotation, 
                                  distance_to_promoters=distance_to_promoters, 
                                  remove_blacklist_regions=TRUE, 
                                  ENCODE_blacklist=ENCODE_blacklist)

#SAVE counts
if(N==1000000){
  N=length(enhancers)
}

regions=enhancers[1:N]

system.time(enhancer_profiles_parallel<-extract_profiles_parallel(bam_folder=bam_folder, 
                                      regions=regions, directionality=FALSE, 
                                      directions=strand(regions), window=window, bin_size=bin_size))


setwd(bam_folder)
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
for(i in 1:length(enhancer_profiles_parallel)){
  profiles_temp[[mark_name[i]]]=enhancer_profiles_parallel[[i]][[1]]
  counts[[mark_name[i]]]=enhancer_profiles_parallel[[i]][[2]]
}

profiles$profiles=profiles_temp
profiles$counts=counts



normalized_profiles=normalize(profiles, round_logic, negatives_to_zero)

#normalize wrt to other cell line if applicaple

if(normalizeBool==TRUE){
  otherCellLine_normalized_profiles=normalize_Between_Cell_Lines(profiles=normalized_profiles, counts=counts, 
                                                      countsOtherCellLine=countsOtherCellLine, mods=names(normalized_profiles), 
                                                      round=FALSE, negatives_to_zero=FALSE)
  rm(counts, profiles_temp, enhancer_profiles_parallel)
  save(profiles, otherCellLine_normalized_profiles, normalized_profiles, regions, file=paste(path,cell_line,"/data_R/",NormCellLine,"_normalized_",N,"_enhancers_bin_",bin_size,"_window_",window,".RData",sep=""))


}else{
  #save counts and enhancer profiles
  rm(counts, profiles_temp, enhancer_profiles_parallel)
  save(profiles, normalized_profiles, regions, file=paste(path,cell_line,"/data_R/",N,"_enhancers_bin_",bin_size,"_window_",window,".RData",sep=""))
  
}



