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


chrSizes = seqlengths(Hsapiens)

option_list = list(
  make_option(c("-window", "--window"), type="integer", default=5000, 
              help="window size [default=%default]", metavar="integer"),
  make_option(c("-binSize", "--binSize"), type="integer", default=100, help="bin size (resolution) [default= %default]", metavar="integer"),
  make_option(c("-type", "--type"), type="character", default="maxscore", 
              help=" [default= %default]", metavar="character"),
  make_option(c("-cellLine", "--cellLine"), type="character", default="K562", 
              help="cell line [default= %default]", metavar="character"),
  make_option(c("-distanceMeasure", "--distanceMeasure"), type="character", default="ML", 
              help="distanceMeasure [default= %default]", metavar="character"),
  make_option(c("-e", "--enhancerSeparation"), type="integer", default=2000, 
              help="minimal distance between two enhancers as bp [default= %default]", metavar="integer"),
  make_option(c("-p", "--distToPromoter"), type="integer", default=2000, 
              help="min distance to any TSS as bp [default= %default]", metavar="integer"),
  make_option(c("-overlap", "--overlap"), type="integer", default=100, 
              help="distance between adjacent windows [default= %default]", metavar="integer"),
  make_option(c("-threshold", "--threshold"), type="double", default=0.5, 
              help="threshold [default= %default]", metavar="double"),
  make_option(c("-randomStr", "--randomStr"), type="character", default="pure_random", 
              help="random type [default= %default]", metavar="character"),
  make_option(c("-pathToDir", "--pathToDir"), type="character", default="", 
              help="path to main folder [default= %default]", metavar="character")
); 


opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);


#if (is.null(opt)){
#  print_help(opt_parser)
#  stop("At least one argument must be supplied (input file).n", call.=FALSE)
#}


window=opt$window								
distance_to_promoters=opt$distToPromoter										
bin_size=opt$binSize
type=opt$type
enhancer_separation=opt$enhancerSeparation
cell_line=opt$cellLine
distance_measure=opt$distanceMeasure
overlap=opt$overlap
threshold=opt$threshold
random_str=opt$randomStr
path_to_dir=opt$pathToDir


predictions_path=paste(path_to_dir,"/results/model_promoters_and_random_combined/",sep="")
results_path=paste(predictions_path,cell_line,"/",
                   "/",distance_measure,"/",sep="")



################################Load windows#########################################################################
load(paste(results_path, "whole_genome_data.RData",sep="")) #38.6 GB

keep(enhancer_regions, negative_regions, type,
     enhancer_separation,
     distance_measure,
     cell_line,
     window,
     bin_size,
     threshold,
     overlap,
     distance_to_promoters,
     predictions_path,
     results_path, GRanges_prediction_regions, path_to_dir, sure=TRUE) #633 MB
#keep enhancer_regions,

###################################Load libsvm predictions #####################################################################################

data=read.table(paste(results_path, "/bin_",bin_size,"_whole_genome_test.txt.predict",sep=""),
                header=TRUE, stringsAsFactors=FALSE) 

colnames(data)<-c("labels", "enhancer_score", "non-enhancer score")
#######################################################################

source("code/functions.R")


predictions_full=unlist(GRangesList(GRanges_prediction_regions), use.names=FALSE) #this creates the dublicate chromosome problem if use.names=TRUE
rm(GRanges_prediction_regions)
predictions_full$label=data$labels
predictions_full$enhancer_score=data$enhancer_score
rm(data)

##write all predictions to bedgraph, do not process anyhow

##########################bedgraph file for all predictions#####################################################################3
#these should be centered to 100 bp windows

tmp=resize(predictions_full, width=100, fix="center" )

tmp=tmp[order(tmp)]
tmp<-as.data.frame(tmp)
names(tmp)=c("seqnames", "start", "end", "width", "strand","label","score")
tmp$start=tmp$start-1

filename=paste(results_path, "prediction_scores_PREPRINT.bedGraph",sep="")

write.table(tmp[,c("seqnames", "start", "end", "score")], file=filename, quote=FALSE, sep=" ", row.names=FALSE, col.names=FALSE)



#check this whether the function is correct, and there is no mess with the bin_size and overlap!!!!
#do not use type="multi", does not work

#the histogram of enhancer lengths could be more informative
predictions_full=enhancer_grouping_choose_type(predictions_full=predictions_full, 
                                               filename=paste(results_path,
                                                              "enhancer_lengths_",sep=""), 
                                               type=type, 
                                               threshold=threshold,
                                               enhancer_separation=enhancer_separation, 
                                               overlap=overlap)

#remove those overlapping ENCODE blacklist regions
ENCODE_blacklist=ENCODE_blaclist_regions(path_to_dir)

test=as.matrix(findOverlaps(predictions_full, ENCODE_blacklist))
if(nrow(test)!=0){
  predictions_full=predictions_full[-unique(test[,1])]
}

order_ind=order(predictions_full$enhancer_score, decreasing=TRUE)

predictions_full=predictions_full[order_ind]

save(predictions_full, file=paste(results_path, "PREPRINT_predictions_",type,"_",threshold,".RData",sep=""))



########################################################################################

enhancer_ind=which(predictions_full$enhancer_score>threshold) 
enhancer_predictions_all<-list()
enhancer_predictions_all[[1]]<-predictions_full[enhancer_ind]

nonenhancer_ind=which(predictions_full$enhancer_score<=threshold) #492156
nonenhancer_predictions_all<-list()
nonenhancer_predictions_all[[1]]<-predictions_full[nonenhancer_ind]



ElementMetadata_enhancer_predictions_all<-list()
ElementMetadata_enhancer_predictions_all[[1]]<-elementMetadata(enhancer_predictions_all[[1]])

ElementMetadata_nonenhancer_predictions_all<-list()
ElementMetadata_nonenhancer_predictions_all[[1]]<-elementMetadata(nonenhancer_predictions_all[[1]])

#####################which are close to TSS##############################################################################
#TSS_annotation:
#GFF/GTF File format, start and end are both one-based

path=paste(path_to_dir,"/Data/",sep="")

GR_Gencode_protein_coding_TSS_positive=readRDS(paste(path,"GENCODE_TSS/","GR_Gencode_protein_coding_TSS_positive.RDS",sep=""))


#these need the chromosome length information

human.chromlens = seqlengths(Hsapiens)

seqlengths(GR_Gencode_protein_coding_TSS_positive)<-human.chromlens[seqnames(seqinfo(GR_Gencode_protein_coding_TSS_positive))]


TSS_annotation=GR_Gencode_protein_coding_TSS_positive

ElementMetadata_enhancer_predictions_all[[1]]$TSS= abs( distToTss(enhancer_predictions_all[[1]],TSS_annotation)  ) <= distance_to_promoters

ElementMetadata_nonenhancer_predictions_all[[1]]$TSS= abs(distToTss(nonenhancer_predictions_all[[1]],TSS_annotation)  ) <= distance_to_promoters


#####################Remove all training regions from the predictions #####################################################
mcols(enhancer_regions)<-NULL
training_regions=c(enhancer_regions, negative_regions)


mtch<-try( abs(distToTss(enhancer_predictions_all[[1]], training_regions)) <= window/2, silent=TRUE) #window=1500

if( (class(mtch)!="try-error") & (length(mtch)!=0) ){
  remove_index=which(mtch==TRUE)
  
  if(length(remove_index)!=0){
    enhancer_predictions_all[[1]]=enhancer_predictions_all[[1]][-remove_index]
    
    
    ElementMetadata_enhancer_predictions_all[[1]]=ElementMetadata_enhancer_predictions_all[[1]][-remove_index,]
    
  }
  
}


mtch<-try( abs(distToTss(nonenhancer_predictions_all[[1]], training_regions)) <= window/2, silent=TRUE) #window=1500

if( (class(mtch)!="try-error") & (length(mtch)!=0) ){
  remove_index=which(mtch==TRUE)
  
  if(length(remove_index)!=0){
    nonenhancer_predictions_all[[1]]=nonenhancer_predictions_all[[1]][-remove_index]
    
    
    ElementMetadata_nonenhancer_predictions_all[[1]]=ElementMetadata_nonenhancer_predictions_all[[1]][-remove_index,]
    
  }
  
}





########################## With and without TSS ####################################################################################


#and now we investigate with TSS and without TSS
#enhancers_with_TSS

enhancer_predictions_with_TSS<-list()
enhancer_predictions_with_TSS[[1]]<-enhancer_predictions_all[[1]][which(ElementMetadata_enhancer_predictions_all[[1]]$TSS==TRUE)]

ElementMetadata_enhancer_predictions_with_TSS<-list()
ElementMetadata_enhancer_predictions_with_TSS[[1]]<-ElementMetadata_enhancer_predictions_all[[1]][which(ElementMetadata_enhancer_predictions_all[[1]]$TSS==TRUE),]


#enhancers_without_TSS

enhancer_predictions_without_TSS<-list()
enhancer_predictions_without_TSS[[1]]<-enhancer_predictions_all[[1]][which(ElementMetadata_enhancer_predictions_all[[1]]$TSS==FALSE)]

ElementMetadata_enhancer_predictions_without_TSS<-list()
ElementMetadata_enhancer_predictions_without_TSS[[1]]<-ElementMetadata_enhancer_predictions_all[[1]][which(ElementMetadata_enhancer_predictions_all[[1]]$TSS==FALSE),]





predictions_all<-list()
predictions_all$enhancers<-enhancer_predictions_all





predictions_with_TSS<-list()
predictions_with_TSS$enhancers<-enhancer_predictions_with_TSS



predictions_without_TSS<-list()
predictions_without_TSS$enhancers<-enhancer_predictions_without_TSS


ElementMetadatas_all<-list()
ElementMetadatas_all$enhancers<-ElementMetadata_enhancer_predictions_all


ElementMetadatas_with_TSS<-list()
ElementMetadatas_with_TSS$enhancers<-ElementMetadata_enhancer_predictions_with_TSS


ElementMetadatas_without_TSS<-list()
ElementMetadatas_without_TSS$enhancers<-ElementMetadata_enhancer_predictions_without_TSS


predictions<-list()
predictions$all<-predictions_all
predictions$with_TSS<-predictions_with_TSS
predictions$without_TSS<-predictions_without_TSS

ElementMetadatas<-list()
ElementMetadatas$all<-ElementMetadatas_all
ElementMetadatas$with_TSS<-ElementMetadatas_with_TSS
ElementMetadatas$without_TSS<-ElementMetadatas_without_TSS


rm(list=c("predictions_all",
          "predictions_with_TSS",
          "predictions_without_TSS",
          "ElementMetadatas_all",
          "ElementMetadatas_with_TSS",
          "ElementMetadatas_without_TSS"
          
))


#order the predictions based on enhancer prediction score

for(i in 1:length(predictions)){ #all with_TSS without_TSS
  print(i)
  for(j in 1:length(predictions[[i]])){ # enhancers nonenhancers super_nonenhancers random
    print(j)
    order_ind=order(elementMetadata(predictions[[i]][[j]][[1]])$enhancer_score, decreasing=TRUE)
    for(k in 1:length(predictions[[i]][[j]])){ #full narrow narrower now these are not of equal size!
      print(k)
      predictions[[i]][[j]][[k]]=predictions[[i]][[j]][[k]][order_ind]
      ElementMetadatas[[i]][[j]][[k]]=ElementMetadatas[[i]][[j]][[k]][order_ind, ]
      
    }
    
    
  }
  
  
}

save(predictions, ElementMetadatas, 
     file=paste(results_path, "PREPRINT_predictions_",type,"_",threshold,".RData",sep=""))

for(predictions_TSS in c("all", "without_TSS")){
  
  tmp=resize( predictions[[predictions_TSS]][["enhancers"]][[1]], width=100, fix="center")
  
  tmp=tmp[order(tmp)]
  tmp<-as.data.frame(tmp)
  names(tmp)=c("seqnames", "start", "end", "width", "strand","label","score")
  tmp$start=tmp$start-1
  write.table(tmp[,c("seqnames", "start", "end", "score")], 
              file=paste(results_path,"enhancers_PREPRINT_",
              predictions_TSS,"_",type,"_",threshold,".bedGraph",sep=""), 
              quote=FALSE, sep=" ", row.names=FALSE, col.names=FALSE)
  
}


