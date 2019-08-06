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
make_option(c("-cellLine", "--cellLine"), type="character", default="", 
            help="cell line [default= %default]", metavar="character"),
make_option(c("-normalize", "--normalize"), type="logical", default=FALSE, 
            help="do we normalize wrt data from other cell line [default= %default]", metavar="logical"),
make_option(c("-NormCellLine", "--NormCellLine"), type="character", default="", 
            help="name of the cell line normalized wrt [default= %default]", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

window=opt$window
bin_size=opt$binSize
N=opt$N
path_to_dir=opt$pathToDir
cell_line=opt$cellLine
normalizeBool=opt$normalize
NormCellLine=opt$NormCellLine

print(window)
print(bin_size)
print(N)
print(path_to_dir)
print(cell_line)
print(normalizeBool)
print(NormCellLine)

setwd(path_to_dir)
source("code/functions.R")
source("code/extract_profiles_parallel.R")

path=paste(path_to_dir, "/Data",sep="")

#load the counts
if(normalizeBool==TRUE){

  load( file=paste(path,"/",NormCellLine,"/data_R/",as.character(N),"_enhancers_bin_",as.character(bin_size),"_window_",as.character(window),".RData",sep="")) #profiles, normalized_profiles, regions,
  countsOtherCellLine=profiles$counts #remove "V2" from the names
  names(countsOtherCellLine)=gsub("V2","",  names(countsOtherCellLine))



load( file=paste(path,"/",cell_line,"/data_R/",NormCellLine,"_normalized_",N,"_enhancers_bin_",bin_size,"_window_",window,".RData",sep="")) #profiles, normalized_profiles, regions,
counts=profiles$counts #remove "V2" from the names
#names(countsOtherCellLine)=gsub("V2","",  names(countsOtherCellLine))
}else{
  load( file=paste(path,"/",cell_line,"/data_R/",N,"_enhancers_bin_",bin_size,"_window_",window,".RData",sep="")) #profiles, normalized_profiles, regions,
  counts=profiles$counts #remove "V2" from the names
}

bedgraph_path=paste(path_to_dir,"/Data/",cell_line,"/intervals_data_",bin_size,"/",sep="")

intervals_bed_path=paste(path_to_dir,"/Data/intervals_bed_",bin_size,"/",sep="") #zero-start



chrnames=names(Hsapiens)[1:23]
human.chromlens = seqlengths(Hsapiens)[1:23]


split_ranges=GRanges()

#read data for all chromosomes
for(chr in chrnames ){
  print(chr)
  if(chr=="chr1"){
    unionBedGraph= read.table(stringsAsFactors = FALSE,header=TRUE,file=paste(bedgraph_path,"all_",chr,".bedGraph",sep=""))
  }
  else{
    unionBedGraph= rbind(unionBedGraph, read.table(stringsAsFactors = FALSE,header=TRUE,file=paste(bedgraph_path,"all_",chr,".bedGraph",sep="")))
  }
  
  df <- read.table(paste(intervals_bed_path,chr,".bed",sep=""),
                   header=F,
                   stringsAsFactors=F)
  
  
  
  header <- c('chr','start','end','id','score','strand')
  names(df) <- header[1:length(names(df))]
  
  if('strand' %in% colnames(df)){
    df$strand <- gsub(pattern="[^+-]+", replacement = '*', x = df$strand)
  }
  
  
  
  split_ranges <- c(split_ranges, with(df, GRanges(chr, IRanges(start, end), strand=strand))) #243 MB
  
  print(object_size(unionBedGraph)) #in total 2.55 GB
  
  
}


#normalize the data first
non_histone_ind=c(grep('Dnase', names(unionBedGraph)),grep('Nsom', names(unionBedGraph)))
pol_ind=grep('Pol2',names(unionBedGraph))
pol_input=grep('Input',names(unionBedGraph))
mod_input=grep('Control',names(unionBedGraph))

modnames<-names(unionBedGraph)[-c(1,2,3,non_histone_ind, pol_ind, pol_input, mod_input)]

#normalize histones
ControlSignalNumReads=as.vector(counts["Control"])

for(mod in modnames){
  
  #substract from every row
  SignalNumReads=as.vector(counts[mod])
  
  unionBedGraph[,mod]=unionBedGraph[,mod]-unionBedGraph[,"Control"]*(SignalNumReads[[1]]/ControlSignalNumReads[[1]])
}

if(cell_line=="K562"){
SignalNumReads=as.vector(counts["Pol2RawData"])
ControlSignalNumReads=as.vector(counts["InputV2"])

unionBedGraph[,"Pol2RawData"]=unionBedGraph[,"Pol2RawData"]-unionBedGraph[,"InputV2"]*(SignalNumReads[[1]]/ControlSignalNumReads[[1]])

}
if(cell_line=="GM12878"){
  SignalNumReads=as.vector(counts["Pol2RawData"])
  ControlSignalNumReads=as.vector(counts["Input"])
  
  unionBedGraph[,"Pol2RawData"]=unionBedGraph[,"Pol2RawData"]-unionBedGraph[,"Input"]*(SignalNumReads[[1]]/ControlSignalNumReads[[1]])
}
#remove inputs, starts and end and chr
unionBedGraph=unionBedGraph[,-c(1,2,3, pol_input, mod_input)]

if(normalizeBool==TRUE){
  
  #normalize wrt to K562
  for(mod in names(unionBedGraph)){
    unionBedGraph[,mod]=(countsOtherCellLine[[mod]]/counts[[mod]])*unionBedGraph[,mod]
  }
}

unionBedGraph=as.matrix(unionBedGraph) #3.64 G

#where is the signal

#convert negative values to zero
unionBedGraph_zero=unionBedGraph
unionBedGraph_zero[which(unionBedGraph_zero<0)]=0 #Or should the negative ones be converted to positive



save(unionBedGraph, unionBedGraph_zero, split_ranges,
     file=paste(path_to_dir,"/Data/",cell_line,"/data_R/whole_genome_coverage.RData",sep=""))
