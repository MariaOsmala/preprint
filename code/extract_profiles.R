extract_profiles<-function(bam_folder, regions, directionality=TRUE, directions, window, bin_size){
  
  #check that the regions do not go over the chromosome limits, remove if they do
  
  
  
  
  roi=resize(regions,window,'center')
  start_problem=which(start(ranges(roi))<1)
  end_problem=which(end(ranges(roi))>  as.numeric(seqlengths(roi)[as.vector(chrom(roi))]) )
  if(length( c(start_problem, end_problem) )!=0 ){
    roi=roi[-c(start_problem, end_problem)]
    if(directionality==TRUE){
      directions=directions[[c(start_problem, end_problem)]]
      
    }
    
  }
  
  
  setwd(bam_folder)
  bai_files=dir(pattern=".bai")
  bam_files=dir(pattern=".bam")
  bam_files=bam_files[which(bam_files %in% bai_files ==FALSE)]
  
  number_of_reads<-list()
  
  
  
  profiles_all<-list()
  ptm <- proc.time()
  
  for(bam in bam_files){ #this could be parallelized
    
    
    bam_file=paste(bam_folder, "/", bam,sep="")
    gal <- import(bam_file)
    number_of_reads[strsplit(bam,".bam")[[1]]]<-length(gal)
    
    
    #explode each roi
    
    if(bin_size!=1){
      
      #ptm <- proc.time()
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
     
      
      
    }else{ 
    
      cover<-coverage(gal)
      
      profiles=aggregate(cover, roi, FUN=as.vector) #10000 x 8 matrix
      
      if(directionality==TRUE){
        
        
        
        if(length(which(directions=="-"))>1){    
          profiles[, which(directions=="-")]= apply(profiles[, which(directions=="-")],2,rev)
        }
        else{
          profiles[, which(directions=="-")]= rev(profiles[, which(directions=="-")])
        }
        
        
        
      }
      
    }
    
   
    
    profiles_all[[strsplit(bam,".bam")[[1]]]]<-profiles
  }
  
  time1 <- proc.time()-ptm
  
  returned_list<-list()
  returned_list$profiles<-profiles_all
  returned_list$counts<-number_of_reads
  returned_list
  
  
}
