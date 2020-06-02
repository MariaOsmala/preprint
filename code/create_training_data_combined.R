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
          make_option(c("-k", "--k"), type="integer", default=5, 
            help="k-fold CV [default= %default]", metavar="integer"),
make_option(c("-distanceMeasure", "--distanceMeasure"), type="character", default="", 
            help="ML or Bayes_estimated_priors [default= %default]", metavar="character"),
make_option(c("-cellLine", "--cellLine"), type="character", default="", 
            help="cell line [default= %default]", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);


cell_line=opt$cellLine
window=opt$window								
bin_size=opt$binSize
N=opt$N
distance_measure=opt$distanceMeasure #"Bayes_estimated_priors" # ML, Bayes_estimated_priors
k=opt$k #k-fold CV

path_to_dir=opt$pathToDir

print(cell_line)
print(window)
print(bin_size)
print(N)
print(distance_measure)
print(k)

print(path_to_dir)


setwd(path_to_dir)
source("code/functions.R")
path=path_to_dir



directionality=FALSE #for promoters
round_logic=FALSE
negatives_to_zero=FALSE
summary="mean"


################################Load training data enhancers #############################################################
load(file=paste(path,"/Data/",cell_line,"/data_R/",N,"_enhancers_bin_",bin_size,"_window_",window,".RData",sep="")) #normalized_profiles, profiles, regions

original_window=nrow(normalized_profiles[[1]])*bin_size

enhancer_profiles=normalized_profiles
enhancer_regions=regions


enhancer_summary=profile_averages(enhancer_profiles)




#change the window size
#enhancer_profiles=change_window(enhancer_profiles, original_window, window, bin_size)
#enhancer_summary=change_window(enhancer_summary, original_window, window, bin_size)


##################Load training data promoters###########################################################################
#profiles_directed,profiles_undirected, normalized_profiles_directed,normalized_profiles_undirected, regions
load(file=paste(path,"/Data/",cell_line,"/data_R/",N,"_promoters_bin_",
                bin_size,"_window_",window,".RData",sep=""))

promoter_profiles=normalized_profiles_undirected
#promoter_profiles=change_window(promoter_profiles, original_window, window, bin_size)
promoter_regions=regions
#promoter_directions=regions$strand
#averages, quantiles

##################Load training data random################################################################################


#pure_random
load(file=paste(path,"/Data/",cell_line,"/data_R/pure_random_" ,N,"_bin_",
                bin_size,"_window_",window,".RData",sep="")) #profiles, normalized_profiles, regions, accepted_GRanges,steps
  
pure_random_profiles=normalized_profiles
pure_random_regions=regions


#random_with_signal
load(file=paste(path,"/Data/",cell_line,"/data_R/",N,"_random_with_signal_bin_",
                bin_size,"_window_",window,".RData",sep="")) #profiles, normalized_profiles, regions, accepted_GRanges,steps
  
random_with_signal_profiles=normalized_profiles
random_with_signal_regions=regions

#combine negatives

random_profiles=list()
for(i in 1:length(pure_random_profiles)){
  random_profiles[[i]]=cbind(pure_random_profiles[[i]], random_with_signal_profiles[[i]])
}

random_regions=c(pure_random_regions, random_with_signal_regions)

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


###############Create cross validation groups#################################################

cv_groups=createCrossValidationGroups(N_pos=length(enhancer_regions), N_promoter=length(promoter_regions),N_random=length(random_regions) ,k=k)


################Choose function############################################################

if(distance_measure=="ML"){
  
  fn=compute_ML
  
}else if(distance_measure=="Bayes_estimated_priors"){
  
  fn=compute_Bayes_estimated_prior
  
}else{
  print("give correct distance measure")
}

train_summaries_pos<-list() #training data averages or medians
train_summaries_promoters<-list()
train_summaries_random<-list()
for(i in 1:k){
  
  if(summary=="mean"){
      train_summaries_pos[[i]]<-sapply( enhancer_profiles, function (x) rowMeans(x[, cv_groups$pos[[i]]$train]) ) #window length x 15
      train_summaries_promoters[[i]]<-sapply( promoter_profiles, function (x) rowMeans(x[, cv_groups$neg_promoter[[i]]$train]) )
      train_summaries_random[[i]]<-sapply( random_profiles, function (x) rowMeans(x[, cv_groups$neg_random[[i]]$train]) )
  }
 
  
  
  
  train_data_pos=compute_distance_two_negatives(profile=enhancer_profiles, 
                                                subset=cv_groups$pos[[i]]$train, 
                                                summary_pos=train_summaries_pos[[i]], 
                                                summary_neg_promoters=train_summaries_promoters[[i]],
                                                summary_neg_random=train_summaries_random[[i]],
                                                fn=fn, distance_measure=distance_measure,
                                                learn_alpha_prior=TRUE, priorgammas_pos=NULL, 
                                                priorgammas_neg_promoters=NULL, 
                                                priorgammas_neg_random=NULL)

    
  if( distance_measure=="ML"){
    
    
    train_data_neg=compute_distance_two_negatives(profile=negative_profiles, 
                                                  subset=c( cv_groups$neg_promoter[[i]]$train, 
                                                            length(promoter_regions)+cv_groups$neg_random[[i]]$train), 
                                                  summary_pos=train_summaries_pos[[i]],
                                                  summary_neg_promoters=train_summaries_promoters[[i]],
                                                  summary_neg_random=train_summaries_random[[i]],
                                                  fn=fn, distance_measure=distance_measure)
    
    test_data_pos=compute_distance_two_negatives(profile=enhancer_profiles, subset=cv_groups$pos[[i]]$test, 
                                                 summary_pos=train_summaries_pos[[i]], 
                                                 summary_neg_promoters=train_summaries_promoters[[i]],
                                                 summary_neg_random=train_summaries_random[[i]],
                                                 fn=fn,distance_measure=distance_measure)
    
    test_data_neg=compute_distance_two_negatives(profile=negative_profiles, subset=c( cv_groups$neg_promoter[[i]]$test, 
                                                                                      length(promoter_regions)+cv_groups$neg_random[[i]]$test),
                                                 summary_pos=train_summaries_pos[[i]], 
                                                 summary_neg_promoters=train_summaries_promoters[[i]],
                                                 summary_neg_random=train_summaries_random[[i]],
                                                 fn=fn, distance_measure=distance_measure)
  }
   
  if(distance_measure=="Bayes_estimated_priors"){
    
 
    
    
    train_data_neg=compute_distance_two_negatives(profile=negative_profiles, 
                                                  c( cv_groups$neg_promoter[[i]]$train, length(promoter_regions)+cv_groups$neg_random[[i]]$train), 
                                                  summary_pos=train_summaries_pos[[i]],
                                                  summary_neg_promoters=train_summaries_promoters[[i]],
                                                  summary_neg_random=train_summaries_random[[i]],
                                                  fn=fn, distance_measure=distance_measure, 
                                                  learn_alpha_prior=FALSE, priorgammas_pos=train_data_pos$priorgammas_pos, 
                                                  priorgammas_neg_promoters=train_data_pos$priorgammas_neg_promoters,
                                                  priorgammas_neg_random=train_data_pos$priorgammas_neg_random)
    
    
    test_data_pos=compute_distance_two_negatives(profile=enhancer_profiles, subset=cv_groups$pos[[i]]$test, 
                                                 summary_pos=train_summaries_pos[[i]],  
                                                 summary_neg_promoters=train_summaries_promoters[[i]],
                                                 summary_neg_random=train_summaries_random[[i]],
                                                 fn=fn,distance_measure=distance_measure, 
                                                 learn_alpha_prior=FALSE, priorgammas_pos=train_data_pos$priorgammas_pos, 
                                                 priorgammas_neg_promoters=train_data_pos$priorgammas_neg_promoters,
                                                 priorgammas_neg_random=train_data_pos$priorgammas_neg_random)
    
    test_data_neg=compute_distance_two_negatives(profile=negative_profiles, subset=c( cv_groups$neg_promoter[[i]]$test, 
                                                                                      length(promoter_regions)+cv_groups$neg_random[[i]]$test), 
                                                 summary_pos=train_summaries_pos[[i]],  
                                                 summary_neg_promoters=train_summaries_promoters[[i]],
                                                 summary_neg_random=train_summaries_random[[i]],
                                                 fn=fn, distance_measure=distance_measure, 
                                                 learn_alpha_prior=FALSE, priorgammas_pos=train_data_pos$priorgammas_pos, 
                                                 priorgammas_neg_promoters=train_data_pos$priorgammas_neg_promoters,
                                                 priorgammas_neg_random=train_data_pos$priorgammas_neg_random)
    
  
  }
    
    ##########################Visualization of the data######################################################################## 
  
  common_path=paste(path, "/results/model_promoters_and_random_combined/",cell_line,"/",distance_measure,"/",k,"-fold_CV_",i,"/NSamples_",N ,"_window_",window,"_bin_",bin_size,"_",k,"fold_cv_",i,sep="")

  
  
    
  
      
     
      
      
                                                                 
      
  
   
    
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
    
    save.image(file=paste(common_path,"_training_data.RData" ,sep=""))
    
    #########################################unnormalized data is saved in libsvm format########################################################################################
    append_bool=FALSE
    for(h in 1:(ncol(test.data)) ){
      string=paste(test.labels[h]," ",sep="")
      for(j in 1:nrow(test.data)){
        string<-paste(string,j,":",test.data[j,h]," ",sep="")
      }
      if(append_bool==FALSE){
        write(string, file=paste(common_path,"_test_data.txt",sep=""), append=FALSE)
        append_bool=TRUE
      }else{
        write(string, file=paste(common_path,"_test_data.txt",sep="") , append=TRUE)
        
      }
    }
    append_bool=FALSE
    for(h in 1:(ncol(train.data)) ){
      string=paste(train.labels[h]," ",sep="")
      for(j in 1:nrow(train.data) ){
        string<-paste(string,j,":",train.data[j,h]," ",sep="")
      }
      if(append_bool==FALSE){
        write(string, file=paste(common_path,"_train_data.txt",sep=""), append=FALSE)
        append_bool=TRUE
      }else{
        write(string, file=paste(common_path,"_train_data.txt",sep=""), append=TRUE)
      }
    }
    
    
}
