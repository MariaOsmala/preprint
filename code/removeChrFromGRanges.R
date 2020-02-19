removeChrFromGRanges <- function(peaksGRanges, chr){
  if(length(which(seqnames(peaksGRanges)==chr))!=0){
    peaksGRanges=peaksGRanges[-which(seqnames(peaksGRanges)==chr)]
    
  }
  peaksGRanges
  
}