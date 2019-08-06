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

human.chromlens = seqlengths(Hsapiens)


args <- commandArgs(trailingOnly = TRUE)

option_list = list(
  make_option(c("-b", "--binSize"), type="integer", default=100, help="bin size (resolution) [default= %default]", metavar="integer"),
  make_option(c("-output", "--output"), type="character", default="", 
              help="output folder [default= %default]", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);


bin_size=opt$binSize
output_folder=opt$output

chr_names=c("chr1","chr2","chr3",  "chr4",  "chr5","chr6",  
                      "chr7",  "chr8",  "chr9",  "chr10", "chr11", "chr12", 
                      "chr13", "chr14", "chr15", "chr16", "chr17","chr18", 
                      "chr19",   "chr20", "chr21", "chr22",   "chrX") 



for(chr in chr_names){
    
    chrlength=human.chromlens[chr]
    
    #this many windows
    interval_nro=ceiling((chrlength-bin_size)/bin_size)
    start=rep(0, interval_nro)   
    end=rep(0, interval_nro) 
    long_vector=seq(0, (interval_nro-1), 1)
    start=long_vector*bin_size+1    
    end=long_vector*bin_size+ bin_size
    
    strandinformation=Rle( strand( rep( '+',length(start) ) ))
    regions=GRanges(seqnames = Rle(rep(chr,length(start)), rep(1, length(start)) ), ranges = IRanges(start=start, end=end),  strand = strandinformation )
    
    export(object=regions, format = "bed", 
           con=paste(output_folder,chr,".bed",sep=""))
    
    
}
