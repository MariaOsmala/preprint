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
make_option(c("-distanceMeasure", "--distanceMeasure"), type="character", default="", 
            help="ML or Bayes_estimated_priors [default= %default]", metavar="character"),
make_option(c("-cellLine", "--cellLine"), type="character", default="", 
            help="cell line [default= %default]", metavar="character"),
make_option(c("-randomStr", "--randomStr"), type="character", default="", 
            help="random type [default= %default]", metavar="character"),
make_option(c("-pathToDir", "--pathToDir"), type="character", default="", 
            help="path to main folder [default= %default]", metavar="character"),
make_option(c("-NormCellLine", "--NormCellLine"), type="character", default="", 
            help="name of the cell line normalized wrt [default= %default]", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);


cell_line=opt$cellLine
window=opt$window								
bin_size=opt$binSize
N=opt$N
path_to_dir=opt$pathToDir
distance_measure=opt$distanceMeasure #"Bayes_estimated_priors" # ML, Bayes_estimated_priors
random_str=opt$randomStr
NormCellLine=opt$NormCellLine

print(cell_line)
print(window)
print(bin_size)
print(N)
print(path_to_dir)
print(distance_measure)
print(random_str)
print(NormCellLine)


overlap=100

setwd(path_to_dir)
source("code/functions.R")
directionality=FALSE #for promoters
round_logic=FALSE
negatives_to_zero=FALSE



################################Load training data enhancers #############################################################
load(file=paste(path_to_dir,"/Data/",NormCellLine,"/data_R/",
                N,"_enhancers_bin_",bin_size,"_window_",window,".RData",sep="")) #normalized_profiles, profiles, regions

original_window=nrow(normalized_profiles[[1]])*bin_size

enhancer_profiles=normalized_profiles
enhancer_regions=regions


#change the window size
enhancer_profiles=change_window(enhancer_profiles, original_window, window, bin_size)



##################Load training data promoters###########################################################################
#profiles_directed,profiles_undirected, normalized_profiles_directed,normalized_profiles_undirected, regions
load(file=paste(path_to_dir,"/Data/",NormCellLine,"/data_R/",N,"_promoters_bin_",bin_size,"_window_",window,".RData",sep=""))

promoter_profiles=normalized_profiles_undirected
promoter_profiles=change_window(promoter_profiles, original_window, window, bin_size)
promoter_regions=regions
#promoter_directions=regions$strand
#averages, quantiles

##################Load training data random################################################################################

print("testing")
print(random_str)

if(random_str=="pure_random"){
  load(file=paste(path_to_dir,"/Data/",NormCellLine,"/data_R/",random_str,"_" ,N,"_bin_",bin_size,"_window_",window,".RData",sep="")) #profiles, normalized_profiles, regions, accepted_GRanges,steps
  
}


if(random_str=="random_with_signal"){
  load(file=paste(path_to_dir,"/Data/",NormCellLine,"/data_R/",N,"_random_with_signal_bin_",bin_size,"_window_",window,".RData",sep="")) #profiles, normalized_profiles, regions, accepted_GRanges,steps
  
}


random_profiles=normalized_profiles
random_profiles=change_window(random_profiles, original_window, window, bin_size)
random_regions=regions

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

train_summaries_pos<-NULL #training data averages or medians
train_summaries_promoters<-NULL
train_summaries_random<-NULL


train_summaries_pos<-sapply( enhancer_profiles, function (x) rowMeans(x) ) #window length x 15
train_summaries_promoters<-sapply( negative_profiles, function (x) rowMeans(x) )
train_summaries_random<-sapply( negative_profiles, function (x) rowMeans(x) )


  
  ############################Load the whole genome data##################################################################################
#Is this already normalized for GM12878 wrt K562?  YES IT IS!
load(file=paste(path_to_dir,"/Data/",cell_line,"/data_R/whole_genome_coverage.RData",sep=""))
rm(unionBedGraph)
  #unionBedGraph_zero, split_ranges
  #unionBedGraph_zero[which(unionBedGraph_zero<0)]=0 This has been done already
unionBedGraph_zero=round_correct(unionBedGraph_zero)
human.chromlens = seqlengths(Hsapiens)

#change split_ranges from bed format to GRanges format

start(split_ranges)=start(split_ranges)+1

#########################Create normalized bigWig files#############################################################
if(distance_measure=="ML" && random_str=="pure_random"){ #create only once
    for(name in dimnames(unionBedGraph_zero)[[2]]){
      tmp=split_ranges
      tmp$score=unionBedGraph_zero[,name]
      seqlengths(tmp)<-human.chromlens[seqnames(seqinfo(tmp))]
      export(tmp, paste(path_to_dir,"/Data/",cell_line,"/coverage_bigWig/",name,"_bin_",bin_size,".bigWig",sep=""),"bigWig")
  }
}

###########################Separate data for different chromosomes##################################################
seqlengths(split_ranges)<-human.chromlens[seqnames(seqinfo(split_ranges))]
split_ranges_list <- split(split_ranges, seqnames(split_ranges))

unionBedGraph_zero_list<-list()
for(chr in names(split_ranges_list)){
    ind=which(seqnames(split_ranges)==chr)
    unionBedGraph_zero_list[[chr]]=unionBedGraph_zero[which(seqnames(split_ranges)==chr),]
}

rm(unionBedGraph_zero)  

  ###################################Compute the probabilities############################################################################
train_data_pos=compute_distance_two_negatives(profile=enhancer_profiles, subset=1:ncol(enhancer_profiles[[1]]), 
                                                summary_pos=train_summaries_pos, summary_neg_promoters=train_summaries_promoters,
                                                summary_neg_random=train_summaries_random,
                                                fn=fn, distance_measure=distance_measure,
                                                learn_alpha_prior=TRUE, priorgammas_pos=NULL, priorgammas_neg_promoters=NULL, priorgammas_neg_random=NULL)

################################################################################################################################
#this many windows
prediction_regions<-function(window, overlap, human.chromlens, chr){
    interval_nro=ceiling((human.chromlens[chr]-window)/overlap)
    start=rep(0, interval_nro)   
    end=rep(0, interval_nro) 
    long_vector=seq(0, (interval_nro-1), 1)
    start=long_vector*overlap+1    
    end=long_vector*overlap+ window
    regions_chrom<-IRanges(start=start, end=end)
    
    strandinformation=Rle( strand( rep( '+',length(start) ) ))
    regions=GRanges(seqnames = Rle(rep(chr,length(start)), rep(1, length(start)) ), ranges = IRanges(start=start, end=end),  strand = strandinformation )
    regions
}
################################################################################################################

data_matrix_list=list()
prediction_regions_list=list()  
    
no_cores=detectCores(logical=FALSE)
registerDoParallel(no_cores)

if( distance_measure=="ML"){
    
    
    train_data_neg=compute_distance_two_negatives(profile=negative_profiles, subset=1:ncol(negative_profiles[[1]]), 
                                                  summary_pos=train_summaries_pos,
                                                  summary_neg_promoters=train_summaries_promoters,
                                                  summary_neg_random=train_summaries_random,
                                                  fn=fn, distance_measure=distance_measure)
    
    #  profiles_all<-list()
    
    
    results_all <- foreach( j=1:length(names(split_ranges_list))) %dopar%
    {
      
      regions<-prediction_regions(window, overlap, human.chromlens, names(split_ranges_list)[j])
      data_matrix<-matrix(0, nrow=length(regions), ncol=45)
      for(i in 1:length(regions)){ #can not do this in this way, separate data for different chromosomes
          #print(i)
          data=unionBedGraph_zero_list[[names(split_ranges_list)[j]]][(1+(i-1)*(overlap/bin_size)):((i-1)*(overlap/bin_size) +(window/bin_size) ),]
    
          data_list=do.call(c, apply(data, 2, list))
   
    
          data_matrix[i,]=compute_distance_two_negatives(profile=data_list, subset=1, 
                                                 summary_pos=train_summaries_pos, 
                                                 summary_neg_promoters=train_summaries_promoters,
                                                 summary_neg_random=train_summaries_random,
                                                 fn=fn,distance_measure=distance_measure)$data
    
      }
    
      result=list()
      result$data_matrix=data_matrix
      result$prediction_regions=regions
      return(result)
      
    }
    
    data_matrix_list<-list()
    prediction_regions_list<-list()
  
    for(i in 1:length(results_all)){
      data_matrix_list[[names(split_ranges_list)[i]]]=results_all[[i]][[1]]
      prediction_regions_list[[ names(split_ranges_list)[i] ]]=results_all[[i]][[2]]
    }
    
    
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
    results_all <- foreach( j=1:length(names(split_ranges_list))) %dopar%
    {
    
      
      regions<-prediction_regions(window, overlap, human.chromlens, names(split_ranges_list)[j])
      data_matrix<-matrix(0, nrow=length(regions), ncol=45)
     
      for(i in 1:length(regions)){
        #print(i)
        data=unionBedGraph_zero_list[[ names(split_ranges_list)[j] ]][(1+(i-1)*(overlap/bin_size)):((i-1)*(overlap/bin_size) +(window/bin_size) ),]
      
        data_list=do.call(c, apply(data, 2, list))
      
      
        data_matrix[i,]=compute_distance_two_negatives(profile=data_list, subset=1, 
                                                     summary_pos=train_summaries_pos, 
                                                     summary_neg_promoters=train_summaries_promoters,
                                                     summary_neg_random=train_summaries_random,
                                                     fn=fn,distance_measure=distance_measure, 
                                                     learn_alpha_prior=FALSE, priorgammas_pos=train_data_pos$priorgammas_pos, 
                                                     priorgammas_neg_promoters=train_data_pos$priorgammas_neg_promoters,
                                                     priorgammas_neg_random=train_data_pos$priorgammas_neg_random)$data
        
        
      
      }
      result=list()
      result$data_matrix=data_matrix
      result$prediction_regions=regions
      return(result)
    }
   
    data_matrix_list<-list()
    prediction_regions_list<-list()
    
    for(i in 1:length(results_all)){
      data_matrix_list[[names(split_ranges_list)[i]]]=results_all[[i]][[1]]
      prediction_regions_list[[ names(split_ranges_list)[i] ]]=results_all[[i]][[2]]
    }
    
  
}
    
###########################Concatenate the data matrix list and regions list ################################################################### 
rm(data_matrix)
data_matrix=do.call(rbind, data_matrix_list)
GRanges_prediction_regions=unlist(prediction_regions_list)

common_path=paste(path_to_dir, "/results/",cell_line,"/",random_str,"/",distance_measure,"/",sep="")




    
#plot_data_boxplots( list(train_data_pos_norm, train_data_neg_norm, test_data_pos_norm, test_data_neg_norm), distance_measure, mod_index=c(1:45), 
#                  modnames=rep(names(enhancer_profiles),3), 
#                  filename=paste(common_path,"_normalized_boxplots.eps",sep=""), yMin=-5, yMax=5 )
    
#######################################################################################################################################################################################    

save.image(paste(common_path, "whole_genome_data.RData",sep=""))    
#########################################unnormalized data is saved in libsvm format########################################################################################
train.data<-as.data.frame( cbind(train_data_pos$data, train_data_neg$data) )
train.labels<-rep(1,(ncol(train_data_pos$data)+ncol(train_data_neg$data)))
train.labels[(ncol(train_data_pos$data)+1):length(train.labels)]<- -1
train.labels=factor(train.labels)


append_bool=FALSE
for(h in 1:(ncol(train.data)) ){
  string=paste(train.labels[h]," ",sep="")
  for(j in 1:nrow(train.data) ){
    string<-paste(string,j,":",train.data[j,h]," ",sep="")
  }
  if(append_bool==FALSE){
    write(string, file=paste(common_path,"bin_",bin_size,"_train_data.txt",sep=""), append=FALSE)
    append_bool=TRUE
  }else{
    write(string, file=paste(common_path,"bin_",bin_size,"_train_data.txt",sep=""), append=TRUE)
  }
}
    
#####################write the whole genome data#################################################################################################
##memory issues, write in parts
parts_nro=floor(nrow(data_matrix)/100)

append_bool=FALSE
for(i in 1:100){
  print(i)
  if(i==100){
    end_ind=nrow(data_matrix)
  }else{
    end_ind=i*parts_nro
  }
  
  
  parts_ind=( (i-1)*(parts_nro)+1 ):end_ind
  ind_length=length( parts_ind ) 
  
  data_matrix_temp=cbind(rep("3", ind_length), 
                         rep(" ", ind_length), 
                         rep("1:", ind_length), data_matrix[parts_ind,1], rep(" ", ind_length),
                         rep("2:", ind_length), data_matrix[parts_ind,2], rep(" ", ind_length),
                         rep("3:", ind_length), data_matrix[parts_ind,3], rep(" ", ind_length),
                         rep("4:", ind_length), data_matrix[parts_ind,4], rep(" ", ind_length),
                         rep("5:", ind_length), data_matrix[parts_ind,5], rep(" ", ind_length),
                         rep("6:", ind_length), data_matrix[parts_ind,6], rep(" ", ind_length),
                         rep("7:", ind_length), data_matrix[parts_ind,7], rep(" ", ind_length),
                         rep("8:", ind_length), data_matrix[parts_ind,8], rep(" ", ind_length),
                         rep("9:", ind_length), data_matrix[parts_ind,9], rep(" ", ind_length),
                         rep("10:", ind_length), data_matrix[parts_ind,10], rep(" ", ind_length),
                         rep("11:", ind_length), data_matrix[parts_ind,11], rep(" ", ind_length),
                         rep("12:", ind_length), data_matrix[parts_ind,12], rep(" ", ind_length),
                         rep("13:", ind_length), data_matrix[parts_ind,13], rep(" ", ind_length),
                         rep("14:", ind_length), data_matrix[parts_ind,14], rep(" ", ind_length),
                         rep("15:", ind_length), data_matrix[parts_ind,15], rep(" ", ind_length),
                         rep("16:", ind_length), data_matrix[parts_ind,16], rep(" ", ind_length),
                         rep("17:", ind_length), data_matrix[parts_ind,17], rep(" ", ind_length),
                         rep("18:", ind_length), data_matrix[parts_ind,18], rep(" ", ind_length),
                         rep("19:", ind_length), data_matrix[parts_ind,19], rep(" ", ind_length),
                         rep("20:", ind_length), data_matrix[parts_ind,20], rep(" ", ind_length),
                         rep("21:", ind_length), data_matrix[parts_ind,21], rep(" ", ind_length),
                         rep("22:", ind_length), data_matrix[parts_ind,22], rep(" ", ind_length),
                         rep("23:", ind_length), data_matrix[parts_ind,23], rep(" ", ind_length),
                         rep("24:", ind_length), data_matrix[parts_ind,24], rep(" ", ind_length),
                         rep("25:", ind_length), data_matrix[parts_ind,25], rep(" ", ind_length),
                         rep("26:", ind_length), data_matrix[parts_ind,26], rep(" ", ind_length),
                         rep("27:", ind_length), data_matrix[parts_ind,27], rep(" ", ind_length),
                         rep("28:", ind_length), data_matrix[parts_ind,28], rep(" ", ind_length),
                         rep("29:", ind_length), data_matrix[parts_ind,29], rep(" ", ind_length),
                         rep("30:", ind_length), data_matrix[parts_ind,30], rep(" ", ind_length),
                         rep("31:", ind_length), data_matrix[parts_ind,31], rep(" ", ind_length),
                         rep("32:", ind_length), data_matrix[parts_ind,32], rep(" ", ind_length),
                         rep("33:", ind_length), data_matrix[parts_ind,33], rep(" ", ind_length),
                         rep("34:", ind_length), data_matrix[parts_ind,34], rep(" ", ind_length),
                         rep("35:", ind_length), data_matrix[parts_ind,35], rep(" ", ind_length),
                         rep("36:", ind_length), data_matrix[parts_ind,36], rep(" ", ind_length),
                         rep("37:", ind_length), data_matrix[parts_ind,37], rep(" ", ind_length),
                         rep("38:", ind_length), data_matrix[parts_ind,38], rep(" ", ind_length),
                         rep("39:", ind_length), data_matrix[parts_ind,39], rep(" ", ind_length),
                         rep("40:", ind_length), data_matrix[parts_ind,40], rep(" ", ind_length),
                         rep("41:", ind_length), data_matrix[parts_ind,41], rep(" ", ind_length),
                         rep("42:", ind_length), data_matrix[parts_ind,42], rep(" ", ind_length),
                         rep("43:", ind_length), data_matrix[parts_ind,43], rep(" ", ind_length),
                         rep("44:", ind_length), data_matrix[parts_ind,44], rep(" ", ind_length),
                         rep("45:", ind_length), data_matrix[parts_ind,45])
  
  write.table(data_matrix_temp, file=paste(common_path,"bin_",bin_size,"_whole_genome_test.txt",sep="") 
              ,quote=FALSE, row.names=FALSE, col.names=FALSE, append=append_bool)
  append_bool=TRUE
  
}




