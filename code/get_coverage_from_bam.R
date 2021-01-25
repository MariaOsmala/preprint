library('GenomicAlignments')

#' compute the number of reads overlapping each bin
#'
#' @param roi 
#' @param gal 
#'
#' @return matrix of length length(roi) x 1, i.e. vector
#' @export
#'
#' @examples

get_coverage_from_bam<-function(roi, gal){
  
  #roi is feature
  
  #gal is reads
  
  #mode can be one of the pre-defined count methods such as 
  #"Union", "IntersectionStrict", or "IntersectionNotEmpty" 
  
  #"Union" : (Default) Reads that overlap any portion of exactly one feature are counted. 
  #Reads that overlap multiple features are discarded. This is the most conservative of the 3 modes.
  
  #"IntersectionStrict" : A read must fall completely "within" the feature to be counted. 
  #If a read overlaps multiple features but falls "within" only one, 
  #the read is counted for that feature. If the read is "within" multiple features, the read is discarded.
  
  #"IntersectionNotEmpty" : A read must fall in a unique disjoint region of a feature to be counted. 
  #When a read overlaps multiple features, the features are partitioned into disjoint intervals. 
  #Regions that are shared between the features are discarded leaving only the unique disjoint regions. 
  #If the read overlaps one of these remaining regions, it is assigned to the feature the unique disjoint region came from.
  
  #inter.feature(Default TRUE) A logical indicating if the counting mode should 
  #be aware of overlapping features. 
  
  #When TRUE (default), reads mapping to multiple features are dropped (i.e., not counted). 
  #When FALSE, these reads are retained and a count is assigned to each feature they map to.
  #There are 6 possible combinations of the mode and inter.feature arguments. 
  
  #mode="Union" & inter.feature=FALSE == countOverlaps(type=any)
  #mode="IntersectionStrict" & inter.feature=FALSE == countOverlaps(type=within)
  
  #‘IntersectionNotEmpty’ does not reduce to a simple countOverlaps because common (shared) 
  #regions of the annotation are removed before counting.
  
  #summarizeOverlaps is from GenomicAlignments package
  
  #output of summarizeOverlaps is RangedSummarizedExperiment from "summarizedExperiment" package, assays is a function from this package
  #gives an integer matrix of lenght(roi) x 1 
  
  assays(summarizeOverlaps(roi, gal, inter.feature=FALSE))$counts
}
