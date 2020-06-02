library(optparse)
library(ShortRead)
library(rtracklayer)
library(bitops)
library(Biostrings)
library(Rsamtools)
library(foreach)
library(doParallel)
library(GetoptLong)
library(circlize)
library(RColorBrewer)
library(BSgenome.Hsapiens.UCSC.hg19)
library(ggplot2)
library(GenomicRanges)
library(pryr)
library(MASS)
library(gridExtra)
library(reshape2)

library(stringr)
library(accelerometry)

library(BSgenome.Hsapiens.UCSC.hg19)

library(gdata)



###################################RFECS and ChromHMM validation###################################################
option_list = list(
  
  make_option(c("-cell", "--cell"), type="character", default="", 
              help="Cell line [default= %default]", metavar="character"),
  make_option(c("-threshold", "--threshold"), type="character", default="05", 
              help="threshold as char [default= %default]", metavar="character"),
  make_option(c("-thresholdNum", "--thresholdNum"), type="double", default=0.5, 
              help="threshold [default= %default]", metavar="double"),
  make_option(c("-randomStr", "--randomStr"), type="character", default="", 
              help="random type [default= %default]", metavar="character"),
  make_option(c("-pathToDir", "--pathToDir"), type="character", default="", 
              help="path to main folder [default= %default]", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

cell_line="K562"
 threshold="05"
 threshold_num=0.5
 threshold_pred="05"
 threshold_pred_num=0.5
 

 path_to_dir="/m/cs/scratch/csb/projects/enhancer_prediction/experiments/RProjects/preprint/"

cell_line=opt$cell
threshold_pred=opt$threshold
threshold_pred_num=opt$thresholdNum
random_str=opt$randomStr
path_to_dir=opt$pathToDir

#RFECS
RFECS_window=100

window=2000
bin_size=100
overlap=100
N=1000

distance_to_promoters=2000
RFECS_results_path=paste(path_to_dir,"/results/RFECS_combined/",cell_line,
                         "/whole_genome_predictions/",
                         cell_line,"_threshold_",threshold_pred,".txt",sep="")


#################Load training data#########################################################
load(paste(path_to_dir,"/results/model_promoters_and_random_combined/",
           cell_line,"/ML/whole_genome_data.RData",sep=""))
keep( threshold_pred,threshold_pred_num, 
     enhancer_regions, negative_regions, RFECS_window,cell_line, path_to_dir,
     window, bin_size, overlap, N, 
     distance_to_promoters, RFECS_results_path, sure=TRUE)
mcols(enhancer_regions)<-NULL
training_regions=c(enhancer_regions, negative_regions)
chrSizes = seqlengths(Hsapiens)
source("code/functions.R")


###################RFECS predictions as GRanges object#####################################
print(RFECS_results_path)
RFECS<-read.table(RFECS_results_path,stringsAsFactors=FALSE)
names(RFECS)<-c("chr", "start", "enhancer_score")
#these are sorted by the location
RFECS_GRanges<-GRanges( seqnames = Rle(RFECS$chr, rep(1, nrow(RFECS) ) ), 
                        ranges = IRanges(start=RFECS$start, end=RFECS$start), 
                        strand = Rle( strand( rep( '+', nrow(RFECS) ) ), rep(1,nrow(RFECS) ) ), 
                        enhancer_score=RFECS$enhancer_score)
order_ind=order(elementMetadata(RFECS_GRanges)$enhancer_score, decreasing=TRUE)
#this makes the regions sorted by the enhancer score, if multiple with the same score, then ordered based on the location
RFECS_GRanges<-RFECS_GRanges[order_ind] 



#####################Remove blacklists####################################

ENCODE_blacklist=ENCODE_blaclist_regions(path_to_dir)

test=as.matrix(findOverlaps(resize(RFECS_GRanges,width=RFECS_window, fix="center", use.names=TRUE), ENCODE_blacklist))
if(nrow(test)!=0){
  RFECS_GRanges=RFECS_GRanges[-unique(test[,1])] 
}


########################################################################
#widen the predictions
#extend RFECS_GRanges to width of window
RFECS_GRanges_RFECS_window=resize(RFECS_GRanges, width=RFECS_window, fix="center", use.names=TRUE) #100
RFECS_GRanges_full=resize(RFECS_GRanges, width=window, fix="center", use.names=TRUE) #2000




######################Check overlap with TSS##########################################################3

GR_Gencode_protein_coding_TSS=readRDS( paste(path_to_dir,"/Data/GENCODE_TSS/","GR_Gencode_protein_coding_TSS.RDS",sep="")) #73271
GR_Gencode_protein_coding_TSS_positive=readRDS(paste(path_to_dir,"/Data/GENCODE_TSS/","GR_Gencode_protein_coding_TSS_positive.RDS",sep=""))

#these need the chromosome length information

human.chromlens = seqlengths(Hsapiens)

seqlengths(GR_Gencode_protein_coding_TSS_positive)<-human.chromlens[seqnames(seqinfo(GR_Gencode_protein_coding_TSS_positive))]
seqlengths(GR_Gencode_protein_coding_TSS)<-human.chromlens[seqnames(seqinfo(GR_Gencode_protein_coding_TSS))]

TSS_annotation=GR_Gencode_protein_coding_TSS_positive

elementMetadata_RFECS_RFECS_window=data.frame(TSS = (abs(distToTss(RFECS_GRanges,TSS_annotation)  ) <= distance_to_promoters)) #38082
elementMetadata_RFECS_full=data.frame(TSS= elementMetadata_RFECS_RFECS_window$TSS)


######################Remove training data##########################################



mtch<-try( abs(distToTss(RFECS_GRanges_full, training_regions)) <= window/2, silent=TRUE) #window=1500


if( (class(mtch)!="try-error") & (length(mtch)!=0) ){
  
  remove_index=which(mtch==TRUE) #1010
  
  RFECS_GRanges_RFECS_window=RFECS_GRanges_RFECS_window[-remove_index]
  RFECS_GRanges_full=RFECS_GRanges_full[-remove_index] #37072
 
  
  elementMetadata_RFECS_RFECS_window<-elementMetadata_RFECS_RFECS_window[-remove_index,1]
  elementMetadata_RFECS_full<-elementMetadata_RFECS_full[-remove_index,1]
 
  
  
}




enhancer_predictions_RFECS<-list()
enhancer_predictions_RFECS[[1]]=RFECS_GRanges_RFECS_window
enhancer_predictions_RFECS[[2]]=RFECS_GRanges_full

#############################enhancers_with_TSS############################################################################################


enhancer_predictions_with_TSS_RFECS<-list()
enhancer_predictions_with_TSS_RFECS[[1]]<-enhancer_predictions_RFECS[[1]][which(elementMetadata_RFECS_RFECS_window==TRUE)]
enhancer_predictions_with_TSS_RFECS[[2]]<-enhancer_predictions_RFECS[[2]][which(elementMetadata_RFECS_full==TRUE)]





#############################enhancers without_TSS######################################################################################3
enhancer_predictions_without_TSS_RFECS<-list()
enhancer_predictions_without_TSS_RFECS[[1]]<-enhancer_predictions_RFECS[[1]][which(elementMetadata_RFECS_RFECS_window==FALSE)]
enhancer_predictions_without_TSS_RFECS[[2]]<-enhancer_predictions_RFECS[[2]][which(elementMetadata_RFECS_full==FALSE)]


##############################nonenhancers#########################################################################################333
predictions_RFECS=list()
predictions_RFECS$all=enhancer_predictions_RFECS[[1]]
predictions_RFECS$with_TSS=enhancer_predictions_with_TSS_RFECS[[1]]
predictions_RFECS$without_TSS=enhancer_predictions_without_TSS_RFECS[[1]]
  
  
save(predictions_RFECS,file=paste(path_to_dir, 
                                  "/results/RFECS_combined/",cell_line,
                                                          "/whole_genome_predictions/",
                                  cell_line,"_threshold_",threshold_pred,".RData",sep=""))
###########################bedfiles################################################

rfecs_bedfiles=paste(path_to_dir,"/results/RFECS_combined/",
                     cell_line,"/bedfiles/",sep="")


rfecs_all_prediction_scores=paste(path_to_dir,"/results/RFECS_combined/",
                                  cell_line,"/whole_genome_predictions/",
                                  cell_line,"_allbins_threshold_",threshold_pred,".txt",sep="") #15

for(predictions_TSS in c("all", "without_TSS")){
  tmp=predictions_RFECS[[predictions_TSS]] #100 bp
  tmp=tmp[order(tmp)]
  tmp<-as.data.frame(tmp)
  names(tmp)=c("seqnames", "start", "end", "width", "strand","score")
  tmp$start=tmp$start-1
  filename=paste(rfecs_bedfiles, 
                 "enhancers_",predictions_TSS,"_RFECS_threshold_",threshold_pred,".bedGraph",sep="")
  write.table(tmp[,c("seqnames", "start", "end", "score")], file=filename, quote=FALSE, sep=" ", row.names=FALSE, col.names=FALSE)
}


#all prediction scores for RFECS

RFECS<-read.table(rfecs_all_prediction_scores,stringsAsFactors=FALSE)
names(RFECS)<-c("chr", "start", "enhancer_score")
#these are sorted by the location

RFECS_GRanges<-GRanges( seqnames = Rle(RFECS$chr, rep(1, nrow(RFECS) ) ), 
                        ranges = IRanges(start=RFECS$start, end=RFECS$start), 
                        strand = Rle( strand( rep( '+', nrow(RFECS) ) ), 
                                      rep(1,nrow(RFECS) ) ), enhancer_score=RFECS$enhancer_score)

RFECS_window=100

RFECS_GRanges_RFECS_window=resize(RFECS_GRanges, width=RFECS_window, fix="center", use.names=TRUE)

#write all enhancers as bedGRaph


tmp=RFECS_GRanges_RFECS_window

ENCODE_blacklist=ENCODE_blaclist_regions(path_to_dir)

test=as.matrix(findOverlaps(tmp, ENCODE_blacklist))
if(nrow(test)!=0){
  tmp=tmp[-unique(test[,1])]
}


tmp<-as.data.frame(tmp)
names(tmp)=c("seqnames", "start", "end", "width", "strand","score")
tmp$start=tmp$start-1
filename=paste(rfecs_bedfiles,"all_prediction_scores_RFECS_threshold_",threshold_pred,".bedGraph",sep="")
write.table(tmp[,c("seqnames", "start", "end", "score")], file=filename, quote=FALSE, sep=" ", row.names=FALSE, col.names=FALSE)
