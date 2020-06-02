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
              help="cell line [default= %default]", metavar="character")
)

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);


bin_size=opt$binSize
window=opt$window
N=opt$N
cell_line=opt$cellLine
path_to_dir=opt$pathToDir

bin_size=100
window=2000
N=1000
cell_line="K562"
path_to_dir="/scratch/cs/csb/projects/enhancer_prediction/experiments/RProjects/preprint"
print(window)				#5000			
print(bin_size)
print(N)                #1000
print(cell_line)

#Training data is only K562 data


path_to_dir=opt$pathToDir


data_path=paste(path_to_dir,"/Data/",cell_line,"/data_R/",sep="")
print(data_path)





output_folder=paste(path_to_dir,"/results/RFECS_combined/",cell_line,"/training/",sep="")
  
if(cell_line=="K562"){
  load(file=paste(data_path,N,"_enhancers_bin_",bin_size,"_window_",window,".RData",sep=""))
}else{
  load(file=paste(data_path,"K562_normalized_",N,"_enhancers_bin_",bin_size,"_window_",window,".RData",sep=""))
}  
    
enhancer_regions=regions #these have *, need to be +
    
load(file=paste(data_path,N,"_promoters_bin_",bin_size,"_window_",window,".RData",sep=""))
    
promoter_regions=regions #have *, need to have +
elementMetadata(promoter_regions)<-NULL
    
#random regions
    
load(file=paste(path_to_dir,"/Data/",cell_line,"/data_R/pure_random_" ,N,"_bin_",bin_size,"_window_",window,".RData",sep="")) #profiles, normalized_profiles, regions, accepted_GRanges,steps
      
pure_random_profiles=normalized_profiles
pure_random_regions=regions

    
  
load(file=paste(path_to_dir,"/Data/",cell_line,"/data_R/",N,"_random_with_signal_bin_",bin_size,"_window_",window,".RData",sep="")) #profiles, normalized_profiles, regions, accepted_GRanges,steps
      
  
random_with_signal_profiles=normalized_profiles
random_with_signal_regions=regions

    
    
random_locations=random_regions=c(pure_random_regions, random_with_signal_regions)#the strand is *, earlier these had + already, BUT WAIT, is the data correctly extracted????
    
#now write different types, these need to be sorted by chromosome location
    
strand(enhancer_regions)="+"
strand(promoter_regions)="+"
strand(random_regions)="+"
    
write.table(as.data.frame(enhancer_regions)[,c(1,2,5)], col.names=FALSE, row.names=FALSE, sep="\t", 
                file=paste(output_folder,"enhancers.txt",sep=""), quote=FALSE)
write.table(as.data.frame(c(promoter_regions,random_locations))[,c(1,2,5)], col.names=FALSE, row.names=FALSE, 
                sep="\t", file=paste(output_folder,"non-enhancers.txt",sep=""), quote=FALSE)
    
    
    
  








