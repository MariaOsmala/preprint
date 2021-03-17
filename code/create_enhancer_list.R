#' Create enhancer list as GRanges object based on p300 binding sites, 
#' DNase I HS peaks which have a distance at least distance_to_promoters
#' to any of TSS of genes given as TSS_annotation. Removes chrM regions. Enhancers are
#' first sorted by signalValue, then by qValue
#' @param p300 p300 peaks as .narrowPeak format, needs to contains qValue and signalValue columns
#' @param DNase_peaks_file DNase I HS peaks as .narrowPeak format
#' @param TSS_annotation TSS as GRanges object, single nucleotide coordinates, the strand needs to be positive
#' @param distance_to_promoters The required distance to any of the genes, most be equal or larger than distance_to_promoters
#' @param remove_blacklist_regions Are the blacklist regions removed from from the list of enhancers, Default: TRUE
#' @param ENCODE_blacklist ENCODE blacklist regions given as GRanges object. If the range of size window centered at enhancers overlaps with the blacklist regions, enhancer is removed
#' @param window Region centered at the enhancer of interest, if even, add 1. Enhancer removed if blacklist regions overlap it.
#'
#' @return Enhancers as GRanges object
#' @export
#'
#' @examples
create_enhancer_list<-function(p300, DNase_peaks_file,  TSS_annotation, 
                               distance_to_promoters, 
                               remove_blacklist_regions=TRUE, 
                               ENCODE_blacklist, window){
  
  
  #p300=p300_peaks_file
  #remove_blacklist_regions=TRUE
  
  #this should be correct
  p300peaks_GRanges=narrowPeak2GRanges(narrowPeakFile=p300, asPeak=TRUE)
  
  #remove chrM locations, should chroms Y and X be removed as well???
  
  p300peaks_GRanges=removeChrFromGRanges(p300peaks_GRanges, "chrM")
  
  #8289 of 25881 enhancers removed if all transcripts used
  #5587 of 25881 enhancers removed if protein coding TSS used

  print("Length of p300 file:")
  print(length(p300peaks_GRanges))
  
  strand(p300peaks_GRanges)="+"

  print("Some distances:")
  print(head(distToTss(p300peaks_GRanges,TSS_annotation), 20))
  close_to_TSS=which(abs(distToTss(p300peaks_GRanges,TSS_annotation) ) < distance_to_promoters) 
  
  print("close_to_TSS:")
  print(length(close_to_TSS))
  
  
  GR_Enhancers=p300peaks_GRanges[-close_to_TSS]

  
  
  #DNase_overlap==TRUE    
  
  DNase_peaks_GRanges=narrowPeak2GRanges(DNase_peaks_file, FALSE)
  
  
  
  
  #which p300 peaks overlap DNase peaks, or should we require some distance between these???     
  test=as.matrix(findOverlaps(GR_Enhancers, DNase_peaks_GRanges)) #how many of p300 overlap? 15496 of 20294 113025(Dnase-seq)
  
  
  GR_Enhancers=GR_Enhancers[test[,1]]
  GR_Enhancers=GR_Enhancers[order(GR_Enhancers$signalValue, decreasing=TRUE)] 
  GR_Enhancers=GR_Enhancers[order(GR_Enhancers$qValue, decreasing=TRUE)] #there are many enhancers with the same q-value, but different signalValue
  
  #if window%%2==0 add 1 to window, extend the enhancer window/2 downstream and window/2 upsteam
  #if window%%2==1, extend the enhancer window (window-1)/2 downstream and upstream
  regions_wider=resize(GR_Enhancers, width=window + (window+1)%%2 , fix="center" )
  
  ##################Remove first ENCODE_blacklists from the data##############################################
  
  
  #remove blacklist regions
  test=as.matrix(findOverlaps(regions_wider, ENCODE_blacklist)) #query, subject
  if(nrow(test)!=0){
    GR_Enhancers=GR_Enhancers[-unique(test[,1])]
  }
  
  GR_Enhancers
  
}
