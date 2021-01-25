ENCODE_blaclist_regions<-function(path_to_dir){
  
  # regions in the human genome that have anomalous, unstructured, high signal/read counts in NGS experiments independent of cell line and type of experiment
  
  #empirical from many ENCODE cell lines (not RNA-seq)
  
  #ultra-high signal artifact regions, 226
  #load blaclist regions
  path=paste0(path_to_dir,"/blacklists/")
  
  DAC<-read.table(paste(path,"wgEncodeDacMapabilityConsensusExcludable.bed.gz",sep=""), stringsAsFactors=FALSE) #411
  names(DAC)=c("chr", "start", "end", "name", "score", "strand")
  #[1] "High_Mappability_island" "Satellite_repeat"       
  #[3] "BSR/Beta"                "Low_mappability_island" 
  #[5] "(CATTC)n"                "LSU-rRNA_Hsa"           
  #[7] "centromeric_repeat"      "ALR/Alpha"              
  #[9] "SSU-rRNA_Hsa"            "telomeric_repeat"       
  #[11] "snRNA"                   "TAR1"                   
  #[13] "ACRO1"                   "chrM"                   
  #[15] "(GAGTG)n"                "(GAATG)n"     
  
  Duke<-read.table(paste(path,"wgEncodeDukeMapabilityRegionsExcludable.bed.gz",sep=""), stringsAsFactors=FALSE) #1649
  names(Duke)=c("chr", "start", "end", "name", "score", "strand")
  #[1] "TAR1"         "chrM"         "(GAATG)n"     "BSR/Beta"     "LSU-rRNA_Hsa"
  #[6] "(CATTC)n"     "(GAGTG)n"     "SSU-rRNA_Hsa" "ALR/Alpha"    "ACRO1"       
  #[11] "HSATII"      
  
  strandinformation=Rle( strand( rep( '+',nrow(DAC) ) ))
  DAC_GRanges=GRanges(seqnames = Rle(DAC$chr, rep(1, length(DAC$chr)) ), ranges = IRanges(start=DAC$start+1, end=DAC$end),  strand = strandinformation )
  
  strandinformation=Rle( strand( rep( '+',nrow(Duke) ) ))
  Duke_GRanges=GRanges(seqnames = Rle(Duke$chr, rep(1, length(Duke$chr)) ), ranges = IRanges(start=Duke$start+1, end=Duke$end),  strand = strandinformation )
  
  pathological=union(DAC_GRanges, Duke_GRanges) #1378
  pathological
  
}
