#normalize wrt to input
normalize <- function(profiles, round_logic=FALSE, negatives_to_zero=FALSE){
  
  #profiles is obtained from extract_profiles
  
  counts=profiles$counts
  histone_names=names(profiles$profiles)
  pol_ind=c(grep("Input",histone_names), grep("Pol2", histone_names))
  pol_names=histone_names[pol_ind]
  histone_names=histone_names[-pol_ind]
  other_ind=c(grep("Dnase",histone_names), grep("Nsome", histone_names))
  DNAse_Nucleosome_names=histone_names[other_ind]
  histone_names=histone_names[-other_ind]
  
  DNAse_Nucleosome=profiles$profiles[DNAse_Nucleosome_names]
  histones=profiles$profiles[histone_names]
  polymerase=profiles$profiles[pol_names]
  
  
  normalized_histones<-normalize_Input_and_substract_Input(profiles=histones, 
                                                           readCounts=counts, 
                                                           input_name=histone_names[1], 
                                                           round=round_logic, negatives_to_zero=negatives_to_zero)
  
  normalized_polymerase<-normalize_Input_and_substract_Input(profiles=polymerase, 
                                                             readCounts=counts, 
                                                             input_name=pol_names[1], 
                                                             round=round_logic, negatives_to_zero=negatives_to_zero)
  
  profiles=c(normalized_histones, normalized_polymerase, DNAse_Nucleosome)
  
  
  profiles
  
  
}


normalize_Input_and_substract_Input<-function(profiles, readCounts, input_name, round=TRUE, negatives_to_zero=TRUE){
  #round (Signal − ControlSignal *(SignalNumReads/ControlSignalNumReads))
  
  #this is done for each modification separately, i.e. the input is normalized wrt to each mod at a time
  #also this is done for each chromosome separately
  
  #this function always substract the input
  
  input_data<-profiles[[input_name]]
  
  profiles<-profiles[-which(names(profiles)==input_name)]
  
  for(mod in names(profiles)){
    
    #substract from every row
    SignalNumReads=as.numeric(readCounts[mod])
    ControlSignalNumReads=as.numeric(readCounts[input_name])
    
    #we have just a single vector
    profiles[[mod]]<-profiles[[mod]] - input_data*(SignalNumReads/ControlSignalNumReads)  
    
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