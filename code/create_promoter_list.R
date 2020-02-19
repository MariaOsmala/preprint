#' Create promoter list as GRanges object based on TSS and DNase I HS peaks 
#' Removes chrM regions. 
#' first sorted by pValue then by signalValue, then by distance to DNAse I HS peak summit
#' Consider only those TSS that do not have any other TSS near-by
#'
#' @param TSS_annotation TSS as GRanges object, single nucleotide coordinates
#' @param TSS_annotation_positive TSS as GRanges object, single nucleotide coordinates, the strand needs to be positive
#' @param DNase_peaks_file DNase I HS peaks as .narrowPeak format
#' @param remove_blacklist_regions Are the blacklist regions removed from from the list of promoters, Default: TRUE
#' @param ENCODE_blacklist ENCODE blacklist regions given as GRanges object
#' @param window Region centered at the promoter of interest, if even, add 1
#'
#' @return TSS overlapping DNase-seq peaks, with DNase signalValues, pValues and distance of TSS to DNase peak summit
#' @export
#'
#' @examples
create_promoter_list<-function(TSS_annotation, TSS_annotation_positive, DNase_peaks_file,
                               between_TSS_distance, 
                               remove_blacklist_regions=TRUE, ENCODE_blacklist, window){
  
  #consider only those TSS that do not have any other TSS near-by
  
  
  
  GR_DNase=narrowPeak2GRanges(narrowPeakFile = DNase_peaks_file, asPeak=FALSE)
  GR_DNase_peak=narrowPeak2GRanges(narrowPeakFile = DNase_peaks_file, asPeak=TRUE)
  strand(GR_DNase_peak)="+"
  strand(GR_DNase)="+"
  
  testPeak=as.data.frame( distanceToNearest(GR_DNase_peak,TSS_annotation_positive))
  
  test=as.data.frame( distanceToNearest(GR_DNase,TSS_annotation_positive))
  
  #sort by the distance
  order_ind=order(test$distance, decreasing=FALSE)
  test=test[order_ind,]
  testPeak=testPeak[order_ind,]
  
  TSS_with_DNase=TSS_annotation_positive[test$subjectHits]
  TSS_with_DNase$DNase_signalValue=elementMetadata(GR_DNase[test$queryHits])$signalValue
  TSS_with_DNase$DNase_pValue=elementMetadata(GR_DNase[test$queryHits])$pValue
  TSS_with_DNase$distance=test$distance
  TSS_with_DNase$Peakdistance=testPeak$distance
  TSS_with_DNase$direction=strand(TSS_annotation[test$subjectHits])
  
  test2=as.matrix(findOverlaps(TSS_with_DNase, GR_DNase))
  
  TSS_with_DNase=TSS_with_DNase[unique(test2[,1])] #60528
  
  i=1
  while(i<=length(TSS_with_DNase)){
    test3=distance(TSS_with_DNase[i], TSS_with_DNase) #distances to the same chromosome
    too_close=which(test3<=between_TSS_distance & is.na(test3)==FALSE) 
    
    if(length(too_close)> 1){
      TSS_with_DNase=TSS_with_DNase[-too_close] #remove everything, also i
      
      
    }
    else{
      i=i+1
    }
    
  }
  
  
  #Results about 3000 TSS
  
  
  
  TSS_with_DNase=TSS_with_DNase[which(TSS_with_DNase$distance==0)]
  
  TSS_with_DNase=TSS_with_DNase[order(TSS_with_DNase$DNase_pValue, decreasing=TRUE)]
  TSS_with_DNase=TSS_with_DNase[order(TSS_with_DNase$DNase_signalValue, decreasing=TRUE)]
  TSS_with_DNase=TSS_with_DNase[order(TSS_with_DNase$Peakdistance, decreasing=FALSE)]
  
  
  regions_wider=resize(TSS_with_DNase, width=window + (window+1)%%2 , fix="center" )
  
  #what about the blacklist regions, are these zero or 1-based
  #blacklist regions are 1-based
  if(remove_blacklist_regions==TRUE){
    test=as.matrix(findOverlaps(TSS_with_DNase, ENCODE_blacklist))
    if(nrow(test)!=0){
      TSS_with_DNase=TSS_with_DNase[-unique(test[,1])]
      
    }
    
  }
  
  
  TSS_with_DNase
  
  
}
