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
              help="path to main folder [default= %default]", metavar="character")
 
 
 
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);


window=opt$window					#5000			
bin_size=opt$binSize      #100
N=opt$N                 #1000

path_to_dir=opt$pathToDir
cell_line=opt$cellLine


library(BSgenome.Hsapiens.UCSC.hg19)
chrSizes = seqlengths(Hsapiens)
window=2000
bin_size=100
N=1000
path_to_dir="/scratch/cs/csb/projects/enhancer_prediction/experiments/RProjects/preprint"

random_str="combined"

path=paste(path_to_dir,"/Data/GM12878/",sep="")
load(file=paste(path,"data_R/K562_normalized_",N,"_enhancers_bin_",bin_size,"_window_",window,".RData",sep="")) #normalized_profiles, profiles, regions


enhancer_regions_test=regions

##################Load test data promoters###########################################################################
#profiles_directed,profiles_undirected, normalized_profiles_directed,normalized_profiles_undirected, regions
load(file=paste(path,"/data_R/K562_normalized_",N,"_promoters_bin_",bin_size,"_window_",window,".RData",sep=""))

promoter_regions_test=regions

##################Load test data random################################################################################


load(file=paste(path,"/data_R/K562_normalized_pure_random_" ,
                N,"_bin_",bin_size,"_window_",window,".RData",sep="")) #profiles, normalized_profiles, regions, accepted_GRanges,steps
  

pure_random_profiles=normalized_profiles
pure_random_regions=regions



load(file=paste(path,"/data_R/K562_normalized_",
                N,"_random_with_signal_bin_",bin_size,"_window_",window,".RData",sep="")) #profiles, normalized_profiles, regions, accepted_GRanges,steps
  
random_with_signal_profiles=normalized_profiles
random_with_signal_regions=regions





#combine random
random_profiles=list()
for(i in 1:length(pure_random_profiles)){
  random_profiles[[i]]=cbind(pure_random_profiles[[i]], random_with_signal_profiles[[i]])
}

random_regions_test=c(pure_random_regions, random_with_signal_regions)


tmp=promoter_regions_test
elementMetadata(tmp)=NULL
negative_regions_test=c(tmp, random_regions_test)  
  

keep(enhancer_regions_test, negative_regions_test, random_str,
     
     cell_line,
     N,
     window,
     bin_size,
     path,path_to_dir,sure=TRUE) #633 MB

#################load RFECS predictions###############################################



RFECS_nonenhancers<-read.table(paste(path_to_dir,"/results/RFECS_combined/GM12878/whole_genome_predictions/GM12878_allbins_threshold_05.txt",sep=""),stringsAsFactors=FALSE)
names(RFECS_nonenhancers)<-c("chr", "start", "enhancer_score")
#these are sorted by the location


RFECS_nonenhancers_GRanges_all<-GRanges( seqnames = Rle(RFECS_nonenhancers$chr, rep(1, nrow(RFECS_nonenhancers) ) ), 
                                         ranges = IRanges(start=RFECS_nonenhancers$start, end=RFECS_nonenhancers$start), 
                                         strand = Rle( strand( rep( '+', nrow(RFECS_nonenhancers) ) ), 
                                                       rep(1,nrow(RFECS_nonenhancers) ) ), 
                                         enhancer_score=RFECS_nonenhancers$enhancer_score) #30 milj


#find RFECS predictions that are closest to the test data coordinates
enhancer_predictions_ind=as.data.frame(distanceToNearest(enhancer_regions_test, RFECS_nonenhancers_GRanges_all))

negative_predictions_ind=as.data.frame(distanceToNearest(negative_regions_test, RFECS_nonenhancers_GRanges_all))

predictions=c(RFECS_nonenhancers_GRanges_all[enhancer_predictions_ind$subjectHits]$enhancer_score, RFECS_nonenhancers_GRanges_all[negative_predictions_ind$subjectHits]$enhancer_score)

true_labels=c(rep(1, length(enhancer_regions_test)), rep(-1, length(negative_regions_test)))
#compute AUCs

library("ROCR")

write.table(predictions,file=paste(path_to_dir,"/results/RFECS_combined/GM12878/predictions.txt",sep=""), 
            quote=FALSE,row.names=FALSE, col.names=FALSE)

write.table(true_labels, file=paste(path_to_dir,"/results/RFECS_combined/GM12878/true_labels.txt",sep=""),
            quote=FALSE,row.names=FALSE, col.names=FALSE)

pred.svm<-prediction(predictions, true_labels)

auROC<-performance(pred.svm,'auc')@y.values[[1]] 



write.table(auROC, file=paste(path_to_dir,"/results/RFECS_combined/GM12878/testdata_AUC.txt", sep=""), quote=FALSE,row.names=FALSE, col.names=FALSE)

