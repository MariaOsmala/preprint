quantiles_func<-function(x){
  result<-quantile(x, probs = c(5, 25, 50, 75, 95)/100 ,na.rm = TRUE)
  result
}


profile_averages<-function(enhancers_modifications){
  enhancers_average<-matrix(0,nrow=length(enhancers_modifications), 
                            ncol=nrow(enhancers_modifications[[1]]))
  for(i in 1:length(enhancers_modifications)){
    enhancers_average[i,]=rowMeans(enhancers_modifications[[i]])
    
  }
  
  enhancers_average
}


profile_quantiles<-function(enhancers_modifications){
  #enhancers_modifications is a list
  
  window=nrow(enhancers_modifications[[1]])
  
  quantiles<-list()
  quantiles[[1]]<-matrix(0, nrow=length(enhancers_modifications), ncol=window)
  quantiles[[2]]<-matrix(0, nrow=length(enhancers_modifications), ncol=window)
  quantiles[[3]]<-matrix(0, nrow=length(enhancers_modifications), ncol=window)
  quantiles[[4]]<-matrix(0, nrow=length(enhancers_modifications), ncol=window)
  quantiles[[5]]<-matrix(0, nrow=length(enhancers_modifications), ncol=window)
  
  
  
  for(i in 1:length(enhancers_modifications)){
    
    temp=apply(enhancers_modifications[[i]], 1,quantiles_func)
    quantiles[[1]][i,]<-temp[1,]
    quantiles[[2]][i,]<-temp[2,]
    quantiles[[3]][i,]<-temp[3,] 
    quantiles[[4]][i,]<-temp[4,] 
    quantiles[[5]][i,]<-temp[5,]
    
    
  }
  
  quantiles
}