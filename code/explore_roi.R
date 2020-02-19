explode_roi<-function(roi2, window, bin){
  
  GRanges(seqnames = Rle(as.character(seqnames(roi2[1])), rep(1, length(window/bin)) ), 
          ranges=IRanges(start=seq(start(roi2[1]), end(roi2[1])-bin+1, by=bin), 
                         end=seq(start(roi2[1])+bin-1, end(roi2[1]), bin) ), strand=Rle( strand( rep( '+',window/bin ) )))
}