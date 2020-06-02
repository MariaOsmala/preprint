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
make_option(c("-distanceMeasure", "--distanceMeasure"), type="character", default="", 
            help="ML or Bayes_estimated_priors [default= %default]", metavar="character"),
make_option(c("-cellLine", "--cellLine"), type="character", default="", 
            help="cell line [default= %default]", metavar="character"),
make_option(c("-normalize", "--normalize"), type="logical", default=FALSE, 
            help="do we normalize wrt data from other cell line [default= %default]", metavar="logical"),
make_option(c("-NormCellLine", "--NormCellLine"), type="character", default="", 
            help="name of the cell line normalized wrt [default= %default]", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);


cell_line=opt$cellLine
window=opt$window								
bin_size=opt$binSize
N=opt$N
distance_measure=opt$distanceMeasure #"Bayes_estimated_priors" # ML, Bayes_estimated_priors
path_to_dir=opt$pathToDir
normalizeBool=opt$normalize
NormCellLine=opt$NormCellLine

print(cell_line)
print(window)
print(bin_size)
print(N)
print(distance_measure)

print(path_to_dir)
print(normalizeBool)
print(NormCellLine)


setwd(path_to_dir)
source("code/functions.R")




directionality=FALSE #for promoters
round_logic=FALSE
negatives_to_zero=FALSE


path=path_to_dir
################################Load training data enhancers #############################################################
if(normalizeBool==TRUE){
  load(file=paste(path,"/Data/",NormCellLine,"/data_R/",N,"_enhancers_bin_",bin_size,"_window_",window,".RData",sep="")) #normalized_profiles, profiles, regions

  original_window=nrow(normalized_profiles[[1]])*bin_size

  enhancer_profiles=normalized_profiles
  enhancer_regions=regions


  enhancer_summary=profile_averages(enhancer_profiles)




  #change the window size
  #enhancer_profiles=change_window(enhancer_profiles, original_window, window, bin_size)
  #enhancer_summary=change_window(enhancer_summary, original_window, window, bin_size)


  ##################Load training data promoters###########################################################################
  #profiles_directed,profiles_undirected, normalized_profiles_directed,normalized_profiles_undirected, regions
  load(file=paste(path,"/Data/",NormCellLine,"/data_R/",N,"_promoters_bin_",bin_size,"_window_",window,".RData",sep=""))

  promoter_profiles=normalized_profiles_undirected
  #promoter_profiles=change_window(promoter_profiles, original_window, window, bin_size)
  promoter_regions=regions
  #promoter_directions=regions$strand
  #averages, quantiles

  ##################Load training data random################################################################################


  #pure_random
  load(file=paste(path,"/Data/",NormCellLine,"/data_R/pure_random_" ,N,"_bin_",bin_size,"_window_",window,".RData",sep="")) #profiles, normalized_profiles, regions, accepted_GRanges,steps
  pure_random_profiles=normalized_profiles
  
  pure_random_regions=regions
  
  
  load(file=paste(path,"/Data/",NormCellLine,"/data_R/",N,"_random_with_signal_bin_",bin_size,"_window_",window,".RData",sep="")) #profiles, normalized_profiles, regions, accepted_GRanges,steps
  
  
  random_with_signal_profiles=normalized_profiles
  random_with_signal_regions=regions
  
  #combine negatives
  
  random_profiles=list()
  for(i in 1:length(pure_random_profiles)){
    random_profiles[[i]]=cbind(pure_random_profiles[[i]], random_with_signal_profiles[[i]])
  }
  
  random_regions=c(pure_random_regions, random_with_signal_regions)
  
  
}


##############positive and integer data############################################
for(i in 1:length(enhancer_profiles)){
  enhancer_profiles[[i]][which(enhancer_profiles[[i]]<0)]=0
  enhancer_profiles[[i]]=round_correct(enhancer_profiles[[i]])
  
  promoter_profiles[[i]][which(promoter_profiles[[i]]<0)]=0
  promoter_profiles[[i]]=round_correct(promoter_profiles[[i]])
  
  random_profiles[[i]][which(random_profiles[[i]]<0)]=0
  random_profiles[[i]]=round_correct(random_profiles[[i]])
  
}

####################COmbine negatives###############################################


negative_profiles=mapply(cbind, promoter_profiles, random_profiles, SIMPLIFY=FALSE)
tmp=promoter_regions
elementMetadata(tmp)=NULL
negative_regions=c(tmp, random_regions)


################Choose function############################################################

if(distance_measure=="ML"){
  
  fn=compute_ML
  
}else if(distance_measure=="Bayes_estimated_priors"){
  
  fn=compute_Bayes_estimated_prior
  
}else{
  print("give correct distance measure")
}



  

train_summaries_pos<-sapply( enhancer_profiles, function (x) rowMeans(x) ) #window length x 15
train_summaries_promoters<-sapply( promoter_profiles, function (x) rowMeans(x) )
train_summaries_random<-sapply( random_profiles, function (x) rowMeans(x) )




################################load test data###############################################################
load(file=paste(path,"/Data/",cell_line,"/data_R/",NormCellLine,"_normalized_",N,"_enhancers_bin_",bin_size,"_window_",window,".RData",sep="")) #normalized_profiles, profiles, regions



enhancer_profiles_test=otherCellLine_normalized_profiles
enhancer_regions_test=regions




##################Load test data promoters###########################################################################
#profiles_directed,profiles_undirected, normalized_profiles_directed,normalized_profiles_undirected, regions
load(file=paste(path,"/Data/",cell_line,"/data_R/",NormCellLine,"_normalized_",N,"_promoters_bin_",bin_size,"_window_",window,".RData",sep=""))

promoter_profiles_test=otherCellLine_normalized_profiles_undirected

promoter_regions_test=regions

##################Load test data random################################################################################

load(file=paste(path,"/Data/",cell_line,"/data_R/",NormCellLine,"_normalized_pure_random_" ,N,"_bin_",bin_size,"_window_",window,".RData",sep="")) #profiles, normalized_profiles, regions, accepted_GRanges,steps
pure_random_profiles_test=otherCellLine_normalized_profiles
pure_random_regions_test=regions

load(file=paste(path,"/Data/",cell_line,"/data_R/",NormCellLine,"_normalized_",N,"_random_with_signal_bin_",bin_size,"_window_",window,".RData",sep="")) #profiles, normalized_profiles, regions, accepted_GRanges,steps
random_with_signal_profiles_test=otherCellLine_normalized_profiles
random_with_signal_regions_test=regions

#combine negatives

random_profiles_test=list()
for(i in 1:length(pure_random_profiles_test)){
  random_profiles_test[[i]]=cbind(pure_random_profiles_test[[i]], 
                                  random_with_signal_profiles_test[[i]])
}

random_regions_test=c(pure_random_regions_test, random_with_signal_regions_test)


##############positive and integer data############################################
for(i in 1:length(enhancer_profiles_test)){
  enhancer_profiles_test[[i]][which(enhancer_profiles_test[[i]]<0)]=0
  enhancer_profiles_test[[i]]=round_correct(enhancer_profiles_test[[i]])
  
  promoter_profiles_test[[i]][which(promoter_profiles_test[[i]]<0)]=0
  promoter_profiles_test[[i]]=round_correct(promoter_profiles_test[[i]])
  
  random_profiles_test[[i]][which(random_profiles_test[[i]]<0)]=0
  random_profiles_test[[i]]=round_correct(random_profiles_test[[i]])
  
}


####################COmbine negatives###############################################


negative_profiles_test=mapply(cbind, promoter_profiles_test, random_profiles_test, SIMPLIFY=FALSE)
tmp=promoter_regions_test
elementMetadata(tmp)=NULL
negative_regions_test=c(tmp, random_regions_test)
  
 
  

train_data_pos=compute_distance_two_negatives(profile=enhancer_profiles, 
                                              subset=1:ncol(enhancer_profiles[[1]]), 
                                              summary_pos=train_summaries_pos, 
                                              summary_neg_promoters=train_summaries_promoters,
                                              summary_neg_random=train_summaries_random,
                                              fn=fn, distance_measure=distance_measure,
                                              learn_alpha_prior=TRUE, 
                                              priorgammas_pos=NULL, 
                                              priorgammas_neg_promoters=NULL, 
                                              priorgammas_neg_random=NULL)

    
if( distance_measure=="ML"){
  
  
  
  train_data_neg=compute_distance_two_negatives(profile=negative_profiles, 
                                                subset=1:ncol(negative_profiles[[1]]), 
                                                  summary_pos=train_summaries_pos,
                                                  summary_neg_promoters=train_summaries_promoters,
                                                  summary_neg_random=train_summaries_random,
                                                  fn=fn, distance_measure=distance_measure)
    
  test_data_pos=compute_distance_two_negatives(profile=enhancer_profiles_test, 
                                               subset=1:ncol(enhancer_profiles_test[[1]]), 
                                                 summary_pos=train_summaries_pos, 
                                                 summary_neg_promoters=train_summaries_promoters,
                                                 summary_neg_random=train_summaries_random,
                                                 fn=fn,distance_measure=distance_measure)
    
  test_data_neg=compute_distance_two_negatives(profile=negative_profiles_test, 
                                               subset=1:ncol(negative_profiles_test[[1]]),
                                                 summary_pos=train_summaries_pos, 
                                                 summary_neg_promoters=train_summaries_promoters,
                                                 summary_neg_random=train_summaries_random,
                                                 fn=fn, distance_measure=distance_measure)
  }
   
if(distance_measure=="Bayes_estimated_priors"){
    
   
    
    train_data_neg=compute_distance_two_negatives(profile=negative_profiles, subset=1:ncol(negative_profiles[[1]]), 
                                                  summary_pos=train_summaries_pos,
                                                  summary_neg_promoters=train_summaries_promoters,
                                                  summary_neg_random=train_summaries_random,
                                                  fn=fn, distance_measure=distance_measure, 
                                                  learn_alpha_prior=FALSE, priorgammas_pos=train_data_pos$priorgammas_pos, 
                                                  priorgammas_neg_promoters=train_data_pos$priorgammas_neg_promoters,
                                                  priorgammas_neg_random=train_data_pos$priorgammas_neg_random)
    
    
    test_data_pos=compute_distance_two_negatives(profile=enhancer_profiles_test, subset=1:ncol(enhancer_profiles_test[[1]]), 
                                                 summary_pos=train_summaries_pos,  
                                                 summary_neg_promoters=train_summaries_promoters,
                                                 summary_neg_random=train_summaries_random,
                                                 fn=fn,distance_measure=distance_measure, 
                                                 learn_alpha_prior=FALSE, priorgammas_pos=train_data_pos$priorgammas_pos, 
                                                 priorgammas_neg_promoters=train_data_pos$priorgammas_neg_promoters,
                                                 priorgammas_neg_random=train_data_pos$priorgammas_neg_random)
    
    test_data_neg=compute_distance_two_negatives(profile=negative_profiles_test, subset=1:ncol(negative_profiles_test[[1]]), 
                                                 summary_pos=train_summaries_pos,  
                                                 summary_neg_promoters=train_summaries_promoters,
                                                 summary_neg_random=train_summaries_random,
                                                 fn=fn, distance_measure=distance_measure, 
                                                 learn_alpha_prior=FALSE, priorgammas_pos=train_data_pos$priorgammas_pos, 
                                                 priorgammas_neg_promoters=train_data_pos$priorgammas_neg_promoters,
                                                 priorgammas_neg_random=train_data_pos$priorgammas_neg_random)
    
  
}
    
##########################Visualization of the data######################################################################## COntinue here

common_path=paste(path, "/results/model_promoters_and_random_combined/",cell_line,"/",distance_measure,"/NSamples_",N ,"_window_",window,"_bin_",bin_size,sep="")


                                                
      
  
   
    
######################Normalize####################################3

train_data_mean=rowMeans(cbind(train_data_pos$data, train_data_neg$data))
train_data_sd=apply(cbind(train_data_pos$data, train_data_neg$data),1,sd)
train_data_pos_norm<-(train_data_pos$data-train_data_mean)/train_data_sd
train_data_neg_norm<-(train_data_neg$data-train_data_mean)/train_data_sd

test_data_pos_norm<-(test_data_pos$data-train_data_mean)/train_data_sd
test_data_neg_norm<-(test_data_neg$data-train_data_mean)/train_data_sd


######################################################################################################################################################################################    
#save samples in libsvm format, normalized or unnormalized version

train.data<-as.data.frame( cbind(train_data_pos$data, train_data_neg$data) )
train.labels<-rep(1,(ncol(train_data_pos$data)+ncol(train_data_neg$data)))
train.labels[(ncol(train_data_pos$data)+1):length(train.labels)]<- -1
train.labels=factor(train.labels)

    
test.data<-as.data.frame( cbind(test_data_pos$data, test_data_neg$data) )
test.labels<-rep(1,(ncol(test_data_pos$data)+ncol(test_data_neg$data)))
test.labels[(ncol(test_data_pos$data)+1):length(test.labels)]<- -1
test.labels=factor(test.labels)

save.image(file=paste(common_path, "_training_data.RData" ,sep=""))
    
#########################################unnormalized data is saved in libsvm format########################################################################################
append_bool=FALSE
for(h in 1:(ncol(test.data)) ){
  string=paste(test.labels[h]," ",sep="")
  for(j in 1:nrow(test.data)){
    string<-paste(string,j,":",test.data[j,h]," ",sep="")
  }
  if(append_bool==FALSE){
    write(string, file=paste(common_path, "_test_data.txt",sep=""), append=FALSE)
    append_bool=TRUE
  }else{
    write(string, file=paste(common_path,"_test_data.txt",sep="") , append=TRUE)
    
  }
}
    
    
common_path=paste(path, "/results/model_promoters_and_random_combined/",NormCellLine,"/",distance_measure,"/NSamples_",N ,"_window_",window,"_bin_",bin_size,sep="")


append_bool=FALSE
for(h in 1:(ncol(train.data)) ){
  string=paste(train.labels[h]," ",sep="")
  for(j in 1:nrow(train.data) ){
    string<-paste(string,j,":",train.data[j,h]," ",sep="")
  }
  if(append_bool==FALSE){
    write(string, file=paste(common_path, "_train_data.txt",sep=""), append=FALSE)
    append_bool=TRUE
  }else{
    write(string, file=paste(common_path,"_train_data.txt",sep=""), append=TRUE)
  }
}
    
