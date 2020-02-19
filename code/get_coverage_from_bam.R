get_coverage_from_bam<-function(roi, gal){
  assays(summarizeOverlaps(roi, gal, inter.feature=FALSE))$counts
  
  
}