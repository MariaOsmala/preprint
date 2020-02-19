normalize_Between_Cell_Lines<-function(profiles, counts, countsOtherCellLine, mods, 
                                       round=FALSE, negatives_to_zero=FALSE){
  
  #profiles are already subtracted by the normalized input
  #normalize the reads
  #(number_of_reads_K562/number_of_reads_GM12878)*(aggregate(cover, roi, FUN=as.vector)
  #Input_profiles_Enhancers[[chr]]*(number_of_reads_GM12878/number_of_reads_Input))
  
  for(mod in mods){
    
    profiles[[mod]]=(countsOtherCellLine[[mod]]/counts[[mod]])*profiles[[mod]]
    
    
    if(negatives_to_zero==TRUE){
      profiles[[mod]][which(profiles[[mod]]<0)]=0
      
    }
    
    if(round==TRUE){
      #window/bin x N
      profiles[[mod]]=round_correct(profiles[[mod]])
    }
  }
  profiles
}