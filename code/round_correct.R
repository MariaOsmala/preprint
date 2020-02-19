round_correct<-function(value){
  
  #value can be a single value, vector or matrix
  #positive and negative values need to be handled separately
  
  
  pos_ind=which(value>=0)
  neg_ind=which(value<0)
  
  value[pos_ind]=trunc(value[pos_ind]+0.5)
  value[neg_ind]=trunc(value[neg_ind]-0.5)
  value
  
}