# library(Rsamtools)
# library(snow)
# library(spp)
# library(accelerometry)
# library(Biostrings)
# library(bitops)
library(BSgenome.Hsapiens.UCSC.hg19)
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
  make_option("--threshold", type="integer", default=10, 
              help="Threshold [default= %default]", metavar="integer"),
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
# 
# #if (is.null(opt)){
# #  print_help(opt_parser)
# #  stop("At least one argument must be supplied (input file).n", call.=FALSE)
# #}
# 
# 
# window=opt$window					#5000			
# bin_size=opt$binSize      #100
# N=opt$N                 #1000
# threshold=opt$threshold #5
# path_to_dir=opt$pathToDir
# cell_line=opt$cellLine
# p300_peaks_file=opt$p300File

window=2000
bin_size=100
N=1000
path_to_dir='~/scratch_cs/csb/projects/enhancer_prediction/aaltorse/Data'
cell_line='K562'
p300_peaks_file='~/scratch_cs/csb/projects/enhancer_prediction/aaltorse/Data/K562/raw_data/wgEncodeAwgTfbsSydhK562P300IggrabUniPk.narrowPeak.gz'
threshold=5

print(window)
print(bin_size)
print(N)
print(threshold)
print(path_to_dir)
print(cell_line)
print(p300_peaks_file)

source("code/functions.R")

setwd(path_to_dir)



#load the counts
path=path_to_dir


load( file=paste(path_to_dir,"/",cell_line,"/data_R/whole_genome_coverage.RData",sep="")) #unionBedGraph, unionBedGraph_zero, split_ranges, accepted_GRanges, steps
rm(unionBedGraph)

start(split_ranges)<-start(split_ranges)+1 #strand is +

unionBedGraph_test=unionBedGraph_zero[,-which(colnames(unionBedGraph_zero)=="Nsome")]

#sum
sums=rowSums(unionBedGraph_test)

accepted=which(sums>=threshold)
accepted_GRanges<-split_ranges[accepted]

human.chromlens = seqlengths(Hsapiens)


sum(width(accepted_GRanges))/sum(as.numeric(human.chromlens)) #3 036 303 846 -> 0.4433607

##################Remove first ENCODE_blacklists from the data##############################################
ENCODE_blacklist=ENCODE_blaclist_regions(path_to_dir)
#blacklist regions are 1-based
test=as.matrix(findOverlaps(accepted_GRanges, ENCODE_blacklist)) #strand is positive for both
if(nrow(test)!=0){
  accepted_GRanges=accepted_GRanges[-unique(test[,1])]
   
}
  

####################Remove p300 BSs and protein coding TSS from the list############################################



p300peaks_GRanges=narrowPeak2GRanges(p300_peaks_file, TRUE)


#remove chrM locations, should chroms Y and X be removed as well???

p300peaks_GRanges=removeChrFromGRanges(p300peaks_GRanges, "chrM")
strand(p300peaks_GRanges)="+"


mtch<-try( abs(distToTss(accepted_GRanges, p300peaks_GRanges)) <= window/2, silent=TRUE) 

if( (class(mtch)!="try-error") & (length(mtch)!=0) ){
  
  remove_index=which(mtch==TRUE)
  
  accepted_GRanges=accepted_GRanges[-remove_index]
}

######################Remove protein coding TSS###################################################################
GR_Gencode_protein_coding_TSS=readRDS( paste(path,"/GENCODE_TSS/","GR_Gencode_protein_coding_TSS.RDS",sep="")) #73271
GR_Gencode_protein_coding_TSS_positive=readRDS(paste(path,"/GENCODE_TSS/","GR_Gencode_protein_coding_TSS_positive.RDS",sep=""))

#these need the chromosome length information


seqlengths(GR_Gencode_protein_coding_TSS_positive)<-human.chromlens[seqnames(seqinfo(GR_Gencode_protein_coding_TSS_positive))]
seqlengths(GR_Gencode_protein_coding_TSS)<-human.chromlens[seqnames(seqinfo(GR_Gencode_protein_coding_TSS))]

TSS_annotation=GR_Gencode_protein_coding_TSS_positive


distance_to_promoters=2000
mtch<-try( abs(distToTss(accepted_GRanges, TSS_annotation)) <= distance_to_promoters, silent=TRUE) #window=1500


if( (class(mtch)!="try-error") & (length(mtch)!=0) ){
  
  remove_index=which(mtch==TRUE)
  
  accepted_GRanges=accepted_GRanges[-remove_index]
  
  
}

accepted_GRanges=reduce(accepted_GRanges)

sum(width(accepted_GRanges))/sum(as.numeric(human.chromlens)) #0.424


########################sample only for random regions where the sum of signal in all bins is min 5################################

windows=c(2000,3000,4000,5000)

for(wi in windows){
print( sum(width( accepted_GRanges[which((width(accepted_GRanges)>=wi)==TRUE)]))/ sum(as.numeric(human.chromlens)))


#  [1] 0.04735523
#  [1] 0.02662455
#  [1] 0.01719851
#  [1] 0.01205615

}


accepted_GRanges_original=accepted_GRanges

accepted_GRanges=accepted_GRanges[which((width(accepted_GRanges)>=window)==TRUE)]

accepted_GRanges_notshrinked=accepted_GRanges

#How large proportion this is of the whole genome
print( sum(width(accepted_GRanges_notshrinked))/sum(as.numeric(human.chromlens))  ) #4,7 %

#######################decrease the region length by window/2 from both ends############################

accepted_GRanges_temp <- resize(accepted_GRanges,  width = width(accepted_GRanges)-window/2+1, fix = 'end')
accepted_GRanges_temp <- resize(accepted_GRanges_temp,  width = width(accepted_GRanges_temp)-window/2, fix = 'start')

accepted_GRanges=accepted_GRanges_temp #strand +
##################################Sample random regions########################################################################
#compute the probability for each range

probs2=width(ranges(accepted_GRanges)) / ( sum(width(ranges(accepted_GRanges))) )
index=sample(1:length(accepted_GRanges), N, replace=TRUE, prob=probs2)



regions=GRanges()


for(i in 1:length(index)){
  
    #first sample a region:  accepted_GRanges[index[i]]
  
    #sample a position within the region
  
    location=sample( start(ranges(accepted_GRanges[index[i]])) : end(ranges(accepted_GRanges [index[i]])) ,1)
    
    
    strandinformation=Rle( strand( rep( '+',1 ) ))
    region=GRanges(seqnames = seqnames(accepted_GRanges[index[i]]), 
                    ranges = IRanges(start=location, end=location),  strand = strandinformation)
    
    
    regions=c(regions, region)

  
}

print(N)

#all have strand +
save( regions, accepted_GRanges, accepted_GRanges_notshrinked, 
      file=paste(path,"/",cell_line,"/data_R/",N,"_randomRegions_with_signal_bin_",bin_size,"_window_",window,".RData",sep=""))


