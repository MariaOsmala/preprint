library('doParallel')

#' Extracts ChIP-seq signal from bam files at specific genomic coordinates
#'
#' @param bam_folder Folder where you have .bam files and .bam.bai files of different modifications
#' @param regions regions as GRanges object
#' @param directionality Are the directions of, for example, promoters, taken into account when generating the final matrixes
#' @param directions factor-Rle of same length as regions, "+", "-" or "*"
#' @param window #window size
#' @param bin_size #the resolution of the signal
#'
#' @return
#' @export
#'
#' @examples
extract_profiles_parallel<-function(bam_folder, regions, directionality=TRUE, directions, window, bin_size){
  
  #check that the regions do not go over the chromosome limits, remove if they do
  
  #If the window is even (parillinen) , window/2 units to the left, 
  #smaller than regions, region will be indexes as window/2+1
  
  #If the window is odd(pariton), region in the middle, 
  #floor(window/2) in the both sites, region will be indexed as ceiling(window/2)
  
  #start(roi[1])+window/2 == regions [1]
  
  roi=resize(regions,width=window,'center')
  start_problem=which(start(ranges(roi))<1)
  end_problem=which(end(ranges(roi))>  as.numeric(seqlengths(roi)[as.vector(seqnames(roi))]) )
  if(length( c(start_problem, end_problem) )!=0 ){
    roi=roi[-c(start_problem, end_problem)]
    if(directionality==TRUE){
      directions=directions[c(start_problem, end_problem)]
      
    }
    
  }
  
  #extract the ChIP-seq feature types
  setwd(bam_folder)
  bai_files=dir(pattern=".bai")
  bam_files=dir(pattern=".bam")
  bam_files=bam_files[which(bam_files %in% bai_files ==FALSE)]
  
  
  no_cores=detectCores(logical=FALSE)
  registerDoParallel(no_cores)
#  profiles_all<-list()
 
  
  profiles_all <- foreach( i=1:length(bam_files)) %dopar%
  { #this could be parallelized
    
    bam=bam_files[i]
    bam_file=paste(bam_folder, "/", bam,sep="")
    gal <- import(bam_file) #GAlignments , import from rtracklayer
    number_of_reads<-length(gal) #sequencing depth
    
    gal_subset=gal[as.matrix(findOverlaps(gal, roi))[,1]]
    
    #if window=2000 and bin_size 1, the enhancer middle index is window/2+1
    
    #if window=2000 and bin_size 100, the enhancer middle (original peak summit) is the start coordinate 
    #of the (window/bin_size)/2+1=11 i.e. start(GRanges_roi[[1]][11])
    
    #The enhancer middle peak would be start(roi) + window/2
    
    #if window=5000 and bin_size 100, 5000 is divided into 100 bp windows, results in 50 bins, 
    #THIS IS WRONG!?: the enhancer middle is between the (5000/100)/2 and (5000/100)/2+1 th feature
    #i.e. between 25th and 26th feature, the dashed line should be drawn between these in heatmaps
    
    
    
    GRanges_roi<-lapply(split(roi,seq(1,length(roi),1)), FUN=explode_roi, window=window, bin=bin_size)
      
    #matrix of window/bin_size x N
    profiles=sapply(GRanges_roi, FUN=get_coverage_from_bam, gal=gal_subset) #gives different regions as rows
      
    if(directionality==TRUE){
        
        
      if(length(which(directions=="-"))>1){    
        profiles[, which(directions=="-")]= apply(profiles[, which(directions=="-")],2,rev)
      }
      else{
        profiles[, which(directions=="-")]= rev(profiles[, which(directions=="-")])
      }
    }
      #time3=proc.time() - ptm
      
      
    #}else{ 
      #alternative way, this can not be binned, you can compute the average afterwards
      #use this when the bin_size==1
      
      
      #ptm <- proc.time()
      #cover<-coverage(gal)
      
      #profiles=aggregate(cover, ranges(roi), FUN=as.vector) #10000 x 8 matrix
      
      #if(directionality==TRUE){
        
        
        
       # if(length(which(directions=="-"))>1){    
      #    profiles[, which(directions=="-")]= apply(profiles[, which(directions=="-")],2,rev)
       # }
      #  else{
       #   profiles[, which(directions=="-")]= rev(profiles[, which(directions=="-")])
      #  }
        
        
        
      #}
      
   # }
    
    #time2=proc.time() - ptm
    
    #user  system elapsed 
    #3.164   0.000   3.238 
    #profiles_all[[strsplit(bam,".bam")[[1]]]]<-profiles
    result=list()
    result$profiles=profiles
    result$counts=number_of_reads
    return(result)
   
  }
  
 
 
  return(profiles_all)
  
}
