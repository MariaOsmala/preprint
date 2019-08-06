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
  make_option(c("-pathToDir", "--pathToDir"), type="character", default="", 
              help="path to main folder [default= %default]", metavar="character"),
  make_option(c("-p300File", "--p300File"), type="character", default="", 
              help="path to p300 peak file [default= %default]", metavar="character"),
  make_option(c("-cellLine", "--cellLine"), type="character", default="", 
              help="cell line [default= %default]", metavar="character"),
  make_option(c("-normalize", "--normalize"), type="logical", default=FALSE, 
              help="do we normalize wrt data from other cell line [default= %default]", metavar="logical"),
  make_option(c("-NormCellLine", "--NormCellLine"), type="character", default="", 
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
N=opt$N
path_to_dir=opt$pathToDir
cell_line=opt$cellLine
p300_peaks_file=opt$p300File
normalizeBool=opt$normalize
NormCellLine=opt$NormCellLine


source("code/functions.R")
source("code/extract_profiles_parallel.R")

if(normalizeBool==TRUE){
  load( file=paste(path_to_dir,"/Data/",NormCellLine,"/data_R/",N,"_enhancers_bin_",bin_size,"_window_",window,".RData",sep="")) #profiles, normalized_profiles, regions,
  countsOtherCellLine=profiles$counts #remove "V2" from the names
  names(countsOtherCellLine)=gsub("V2","",  names(countsOtherCellLine))
}


#############################Sample random regions###################################################
#compute the probability for each chromosome, include only chr1-21 and chrX

human.chromlens = seqlengths(Hsapiens)

allowed_chroms=c("chr1",  "chr2",  "chr3",  "chr4",  "chr5",  "chr6",  "chr7",  "chr8",  "chr9", 
"chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18",
"chr19", "chr20", "chr21", "chr22", "chrX") 

human.chromlens=human.chromlens[allowed_chroms]

chrom_probs=as.numeric(human.chromlens)/sum(as.numeric(human.chromlens))

random_chroms=sample(names(human.chromlens), 10000, replace=TRUE, prob=chrom_probs)

#sample random locations, with distance 2.5 from the ends of chromosomes

allowed_ranges=as.numeric(human.chromlens[random_chroms])-window/2
#######################################################
sample_own<-function(x, n=1, window){
  
  sample( (window/2):x , n)
}
###########################################################
locations=sapply(allowed_ranges, sample_own, n=1, window=window)

strandinformation=Rle( strand( rep( '+',length(locations) ) ))
regions=GRanges(seqnames = random_chroms, 
                ranges = IRanges(start=locations, end=locations),  strand = strandinformation)

##############should be expanded #########################

regions_wider=resize(regions, width=window, fix="center" )

##################Remove first ENCODE_blacklists from the data##############################################
ENCODE_blacklist=ENCODE_blaclist_regions(path_to_dir)
#blacklist regions are 1-based
test=as.matrix(findOverlaps(regions_wider, ENCODE_blacklist))
if(nrow(test)!=0){
  regions=regions[-unique(test[,1])]
  regions_wider=regions_wider[-unique(test[,1])]
   
}
  ####################Remove p300 BSs and protein coding TSS from the list############################################


path=paste(path_to_dir,"/Data/",sep="")

p300peaks_GRanges=narrowPeak2GRanges(p300_peaks_file, TRUE, "qValue")


#remove chrM locations, should chroms Y and X be removed as well???

p300peaks_GRanges=removeChrFromGRanges(p300peaks_GRanges, "chrM")



mtch<-try( abs(distToTss(regions, p300peaks_GRanges)) <= 5000/2, silent=TRUE) 

if( (class(mtch)!="try-error") & (length(mtch)!=0) ){
  
  remove_index=which(mtch==TRUE)
  
  regions=regions[-remove_index]
}

######################Remove protein coding TSS###################################################################
GR_Gencode_protein_coding_TSS=readRDS( paste(path,"GENCODE_TSS/","GR_Gencode_protein_coding_TSS.RDS",sep="")) #73271
GR_Gencode_protein_coding_TSS_positive=readRDS(paste(path,"GENCODE_TSS/","GR_Gencode_protein_coding_TSS_positive.RDS",sep=""))

#these need the chromosome length information

human.chromlens = seqlengths(Hsapiens)

seqlengths(GR_Gencode_protein_coding_TSS_positive)<-human.chromlens[seqnames(seqinfo(GR_Gencode_protein_coding_TSS_positive))]
seqlengths(GR_Gencode_protein_coding_TSS)<-human.chromlens[seqnames(seqinfo(GR_Gencode_protein_coding_TSS))]

TSS_annotation=GR_Gencode_protein_coding_TSS_positive


distance_to_promoters=2000
mtch<-try( abs(distToTss(regions, TSS_annotation)) <= distance_to_promoters, silent=TRUE) #window=1500


if( (class(mtch)!="try-error") & (length(mtch)!=0) ){
  
  remove_index=which(mtch==TRUE)
  
  regions=regions[-remove_index]
  
  
}









################################3extract profiles#####################################################################3
regions=regions[1:N]



figure_path=paste(path_to_dir,"/figures/",sep="")

bam_folder=paste(path,cell_line,"/bam_shifted" ,sep="")
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


#save counts and enhancer profiles
if(normalizeBool==FALSE){
save(profiles, normalized_profiles, regions,  
     file=paste(path,cell_line,"/data_R/pure_random_",N,"_bin_",bin_size,"_window_",window,".RData",sep=""))

for(i in 1:length(profiles$profiles)){
  
  profiles$profiles[[i]]=profiles$profiles[[i]][,1:1000]
  
}

for(i in 1:length(normalized_profiles)){
    normalized_profiles[[i]]=normalized_profiles[[i]][,1:1000]
}

regions=regions[1:1000]

save(profiles, normalized_profiles, regions,  file=paste(path,cell_line,"/data_R/pure_random_1000_bin_",bin_size,"_window_",window,".RData",sep=""))

}else{
  
  otherCellLine_normalized_profiles=normalize_Between_Cell_Lines(profiles=normalized_profiles, counts=counts, 
                                                                 countsOtherCellLine=countsOtherCellLine, mods=names(normalized_profiles), 
                                                                 round=FALSE, negatives_to_zero=FALSE)
  rm(counts, profiles_temp, random_profiles)
  
  
  save(profiles, normalized_profiles, otherCellLine_normalized_profiles,regions, 
       file=paste(path,cell_line,"/data_R/",NormCellLine,"_normalized_","pure_random_",N,"_bin_",bin_size,"_window_",window,".RData",sep=""))
  
  for(i in 1:length(profiles$profiles)){
    
    profiles$profiles[[i]]=profiles$profiles[[i]][,1:1000]
    
  }
  
  for(i in 1:length(normalized_profiles)){
    normalized_profiles[[i]]=normalized_profiles[[i]][,1:1000]
    otherCellLine_normalized_profiles[[i]]=otherCellLine_normalized_profiles[[i]][,1:1000]
  }
  
  regions=regions[1:1000]
  
  save(profiles, normalized_profiles, otherCellLine_normalized_profiles, regions,  
       file=paste(path,cell_line,"/data_R/",NormCellLine,"_normalized_","pure_random_1000_bin_",bin_size,"_window_",window,".RData",sep=""))
  
  
}
