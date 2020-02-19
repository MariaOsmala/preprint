quantiles_func<-function(x){
  result<-quantile(x, probs = c(5, 25, 50, 75, 95)/100 ,na.rm = TRUE)
  result
}


profile_averages<-function(enhancers_modifications){
  enhancers_average<-matrix(0,nrow=length(enhancers_modifications), ncol=nrow(enhancers_modifications[[1]]))
  for(i in 1:length(enhancers_modifications)){
    enhancers_average[i,]=rowMeans(enhancers_modifications[[i]])
    
  }
  
  enhancers_average
}
