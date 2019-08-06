extract_profiles_parallel<-function(bam_folder, regions, directionality=TRUE, directions, window, bin_size){
  
  #check that the regions do not go over the chromosome limits, remove if they do
  
  
  #If the window is even (parillinen) , window/2 units to the left, smaller than regions, region will be indexes as window/2+1
  
  #If the window is odd(pariton), region in the middle, floor(window/2) in the both sites, region will be indexed as ceiling(window/2)
  
  roi=resize(regions,window,'center')
  start_problem=which(start(ranges(roi))<1)
  end_problem=which(end(ranges(roi))>  as.numeric(seqlengths(roi)[as.vector(seqnames(roi))]) )
  if(length( c(start_problem, end_problem) )!=0 ){
    roi=roi[-c(start_problem, end_problem)]
    if(directionality==TRUE){
      directions=directions[c(start_problem, end_problem)]
      
    }
    
  }
  
  
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
    gal <- import(bam_file)
    number_of_reads<-length(gal)
    
    
   
    gal_subset=gal[as.matrix(findOverlaps(gal, roi))[,1]]
    
  
    
    GRanges_roi<-lapply(split(roi,seq(1,length(roi),1)), FUN=explode_roi, window=window, bin=bin_size)
      
    profiles=sapply(GRanges_roi, FUN=get_coverage_from_bam, gal=gal_subset) #gives different regions as rows
      
    if(directionality==TRUE){
        
        
      if(length(which(directions=="-"))>1){    
        profiles[, which(directions=="-")]= apply(profiles[, which(directions=="-")],2,rev)
      }
      else{
        profiles[, which(directions=="-")]= rev(profiles[, which(directions=="-")])
      }
    }
    
    result=list()
    result$profiles=profiles
    result$counts=number_of_reads
    return(result)
   
  }
  
 
 
  return(profiles_all)
  
}
