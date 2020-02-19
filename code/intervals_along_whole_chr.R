intervals_along_whole_chr<-function(chr_length, bin_size){
  bin_number=floor( chr_length/bin_size)
  starts<-1+(seq(1,bin_number,1)-1)*bin_size
  ends=1+( seq(1,bin_number,1) )*bin_size-1
  ends[bin_number]=chr_length
  
  
  
  strandinformation<-rep("+",bin_number)
  GRanges(seqnames = Rle(chr,  bin_number ), ranges = IRanges(start=starts, end = ends), 
          strand = Rle(strand(strandinformation), rep(1,bin_number) ) )
  
}