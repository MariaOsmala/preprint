change_window<-function(profiles, original_window, window, bin_size){
  #enhancer_profiles and promoter_profiles are lists
  
  start=(original_window/bin_size)/2-(window/bin_size)/2+1
  
  end=(original_window/bin_size)/2+(window/bin_size)/2
  
  if(class(profiles)=="list"){ #  window length x samples
    
    #what are the indeces that belong to window
    
    
    test=lapply(profiles,function (x) x[start:end,])
    
  }
  #averages and median are matrixes
  if(class(profiles)=="matrix"){
    test=profiles[, start:end]
    
    
  }
  test
  
  
  
  
}