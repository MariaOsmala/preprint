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

print(window)				#5000			
print(bin_size)
print(N)                #1000
print(cell_line)



path_to_dir=opt$pathToDir


data_path=paste(path_to_dir,"/Data/",cell_line,"/data_R/",sep="")
print(data_path)

random_strs=c("pure_random",  "random_with_signal")

for(random_str in random_strs){

output_folder=paste(path_to_dir,"/results/RFECS/",cell_line,"/",random_str,"/training/",sep="")


load(file=paste(data_path,N,"_enhancers_bin_",bin_size,"_window_",window,".RData",sep=""))



enhancer_regions=regions

load(file=paste(data_path,N,"_promoters_bin_",bin_size,"_window_",window,".RData",sep=""))

promoter_regions=regions$promoters
elementMetadata(promoter_regions)<-NULL
                                                      


#random regions

if(random_str=="pure_random"){
  load(file=paste(path_to_dir,"/Data/",cell_line,"/data_R/",random_str,"_" ,N,"_bin_",bin_size,"_window_",window,".RData",sep="")) #profiles, normalized_profiles, regions, accepted_GRanges,steps
  
}

if(random_str=="random_with_signal"){
  load(file=paste(path_to_dir,"/Data/",cell_line,"/data_R/",N,"_random_with_signal_bin_",bin_size,"_window_",window,".RData",sep="")) #profiles, normalized_profiles, regions, accepted_GRanges,steps
  
}



random_locations=regions

#now write different types, these need to be sorted by chromosome location



write.table(as.data.frame(enhancer_regions)[,c(1,2,5)], col.names=FALSE, row.names=FALSE, sep="\t", 
            file=paste(output_folder,"enhancers.txt",sep=""), quote=FALSE)
write.table(as.data.frame(c(promoter_regions,random_locations))[,c(1,2,5)], col.names=FALSE, row.names=FALSE, 
            sep="\t", file=paste(output_folder,"non-enhancers.txt",sep=""), quote=FALSE)



}



