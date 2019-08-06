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
make_option(c("-cellLine", "--cellLine"), type="character", default="", 
       help="cell line [default= %default]", metavar="character")
)

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);


cell_line=opt$cellLine
window=opt$window								
bin_size=opt$binSize
N=opt$N
K=opt$k #k-fold CV
path_to_dir=opt$pathToDir
cell_line=opt$cellLine

distance_measures=c("Bayes_estimated_priors/",
                    "ML/") 


distance_names=c("Bayes estimated priors",
                 "ML") 

random_types=c( "pure_random","random_with_signal")



auROC_matrix<-matrix(0, nrow=1, ncol=length(distance_measures))

auROC_matrices<-list()

colnames(auROC_matrix)=distance_measures

for(random_type in random_types){

  path=paste(path_to_dir,"/results/",cell_line,"/",random_type,"/",sep="")
  
  
  
  for(distance_measure in distance_measures){ 
    rm(predictions)
    rm(true_labels)
    predictions<-vector()
    true_labels<-vector()
    
    for(k in 1:K){
      print(k)
      
      file=paste(path, distance_measure, "5-fold_CV_",k,
                 "/NSamples_",N,"_window_",window,"_bin_",bin_size,"_5fold_cv_",k,"_test_data.txt.predict",sep="")
      
      predictions=c(predictions,read.table(file, header=TRUE)[,2])
      file=paste(path, distance_measure, "5-fold_CV_", k,
                 "/NSamples_",N,"_window_",window,"_bin_",bin_size,"_5fold_cv_",k,"_test_data.txt",sep="")
      true_labels=c(true_labels, read.table(file)[,1])
    }
    
    pred.svm<-prediction(predictions, true_labels)
    perf.svm <- performance(pred.svm, 'tpr', 'fpr')
    auROC_matrix[1,distance_measure]<-performance(pred.svm,'auc')@y.values[[1]] 
    
    
  }
  
  
  auROC_matrices[[random_type]]=auROC_matrix
  print(auROC_matrix)
  save(auROC_matrix, file=paste(path, "auROCs_",bin_size,".RData",sep=""))
  
}
