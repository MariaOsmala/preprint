#' Convert .narrowPeak file format to GRanges object, the regions in GRanges are ordered based on elementMetadata_significance
#'
#' @param narrowPeakFile Filename possibly including the whole path as text, the narrowPeak file needs to have 10 columns
#' @param asPeak Are the point sources called instead of ranges, Default: FALSE
#' 
#'
#' @return
#' @export
#'
#' @examples
narrowPeak2GRanges<-function(narrowPeakFile, asPeak=FALSE){
  
  #This format is used to provide called peaks of signal enrichment based on pooled, normalized (interpreted) data.
  #narrowPeak format == BED6+4 format, zero based
  #chromStart: The first base in a chromosome is numbered 0
  #chromEnd: the first 100 bases of a chromosome are defined aschromStart=0, chromEnd=100, and span the bases numbered 0-99.
  #Indicates how dark the peak will be displayed in the browser (0-1000). If all scores were '0' when the data were submitted to the DCC, 
  #the DCC assigned scores 1-1000 based on signal value. Ideally the average signalValue per base spread is between 100-1000.
  #signalValue - Measurement of overall (usually, average) enrichment for the region.
  #pValue - Measurement of statistical significance (-log10). Use -1 if no pValue is assigned.
  #qValue - Measurement of statistical significance using false discovery rate (-log10). Use -1 if no qValue is assigned.
  
  
  
  extraCols_narrowPeak <- c(signalValue = "numeric", pValue = "numeric",
                            qValue = "numeric", peak = "integer")
  
  peaksGRanges <- rtracklayer::import(narrowPeakFile, format = "BED",
                          extraCols = extraCols_narrowPeak)
  
  library(BSgenome.Hsapiens.UCSC.hg19)
  human.chromlens = GenomeInfoDb::seqlengths(BSgenome.Hsapiens.UCSC.hg19::Hsapiens)
  seqlengths(peaksGRanges)<-human.chromlens[seqnames(seqinfo(peaksGRanges))]
  
  
  
  if(asPeak==TRUE){
    #peak - Point-source called for this peak; 0-based offset from chromStart. Use -1 if no point-source called.
    #the the coordinate is obtained by adding peak to start (this works for GRanges, 
    #single peak coordinate has always separate start and end which is [start + peak -- start+peak+1])
    
    peaksGRanges <- shift(peaksGRanges, peaksGRanges$peak) #shift start and end both by peak
    width(peaksGRanges) <- 1 #resize to start
    elementMetadata(peaksGRanges)[,c("name", "score", "peak")]=NULL
    #peaksGRanges=GRanges(seqnames = Rle(peaks$chrom, rep(1, length(peaks$chrom)) ), 
    #                            ranges = IRanges(start=peaks$chromStart+peaks$peak+1, end=peaks$chromStart+peaks$peak+1),  
    #                            strand = Rle( strand( rep( '+',nrow(peaks) ) )) )
  }else{
    
    
    #peaksGRanges=GRanges(seqnames = Rle(peaks$chrom, rep(1, length(peaks$chrom)) ), 
    #                     ranges = IRanges(start=peaks$chromStart+1, end=peaks$chromEnd),  
    #                     strand = Rle( strand( rep( '+',nrow(peaks) ) )) )
    
  }
  
  peaksGRanges  
}

