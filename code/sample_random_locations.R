sample_random_locations<-function(chr,human.chromlens,widths){
  
  random_location=data.frame(chrom=chr, start=rep(0, length(widths)), end=rep(0, length(widths)))
  
  for(chrname in c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX")){
    #print(chrname)
    chr_ind=which(chr==chrname)
    location<-sample(1:(human.chromlens[chrname]),length(chr_ind), replace=FALSE)
    ind=which(((location-round(widths[chr_ind]/2))<0) || ((location+round(widths[chr_ind]/2))>human.chromlens[chrname]) )
    while(length(ind)!=0){
      location[ind]=sample(1:(human.chromlens[chrname]),length(ind), replace=FALSE)
      ind=which(((location-round(widths[chr_ind]/2))<0) || ((location+round(widths[chr_ind]/2))>human.chromlens[chrname]) )
    }
    
    random_location$start[chr_ind]=location-round(widths[chr_ind]/2)
    random_location$end[chr_ind]=location+(round(widths[chr_ind]/2)-1)
  }
  
  
  random_location
  
}