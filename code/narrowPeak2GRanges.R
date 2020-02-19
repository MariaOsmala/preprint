narrowPeak2GRanges<-function(narrowPeakFile, asPeak=FALSE, elementMetadata_significance="qValue"){
  
  
  #narrowPeak format == BED6+4 format, zero based
  library(BSgenome.Hsapiens.UCSC.hg19)
  human.chromlens = seqlengths(Hsapiens)
  
  peaks=read.table(narrowPeakFile, stringsAsFactors=FALSE)    
  names(peaks)<-c("chrom","chromStart", "chromEnd", "name", "score","strand", "signalValue", "pValue","qValue", "peak")
  
  if(asPeak==TRUE){
    #peak - Point-source called for this peak; 0-based offset from chromStart. Use -1 if no point-source called.
    peaksGRanges=GRanges(seqnames = Rle(peaks$chrom, rep(1, length(peaks$chrom)) ), 
                         ranges = IRanges(start=peaks$chromStart+peaks$peak+1, end=peaks$chromStart+peaks$peak+1),  
                         strand = Rle( strand( rep( '+',nrow(peaks) ) )) )
  }else{
    
    
    peaksGRanges=GRanges(seqnames = Rle(peaks$chrom, rep(1, length(peaks$chrom)) ), 
                         ranges = IRanges(start=peaks$chromStart+1, end=peaks$chromEnd),  
                         strand = Rle( strand( rep( '+',nrow(peaks) ) )) )
    
  }
  
  #for DNase-seq this is pValue, for p300 signalValue
  if(elementMetadata_significance=="qValue"){
    elementMetadata(peaksGRanges) <- data.frame(peaks$qValue, peaks$signalValue)
  }else{
    elementMetadata(peaksGRanges) <- data.frame(peaks$pValue, peaks$signalValue)
  }
  
  seqlengths(peaksGRanges)<-human.chromlens[seqnames(seqinfo(peaksGRanges))]
  
  peaksGRanges  
  
  
  
  
  
  
  
  
  
}
