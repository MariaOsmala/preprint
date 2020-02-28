#' Separates GRanges roi2 into multiple (N=window/bin) smaller ranges of length bin
#'
#' @param roi2 genomic interval as GRanges, needs to be of width window
#' @param window width of the genomic interval
#' @param bin resolution or bin size
#'
#' @return
#' @export
#'
#' @examples
explode_roi<-function(roi2, window, bin){
  
  GRanges(seqnames = Rle(as.character(seqnames(roi2[1])), rep(1, length(window/bin)) ), 
          ranges=IRanges(start=seq(start(roi2[1]), end(roi2[1])-bin+1, by=bin), 
                         end=seq(start(roi2[1])+bin-1, end(roi2[1]), bin) ), 
          strand=Rle( strand( rep( '*',window/bin ) )))
}
