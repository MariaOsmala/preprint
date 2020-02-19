compute_features<-function(profile, subset, summary, fn, distance_measure){
  
  if(distance_measure=="ML"){
    data_tmp<-matrix(0, nrow=length(profile), ncol=length(subset) ) #15 x N
    alphas<-matrix(0, nrow=length(profile), ncol=length(subset) )  #15 x N
    for(p in 1:length(profile)){
      tmp=apply( profile[[p]][,subset], 2, function (x) fn(x,summary[,p]))
      alphas[p,]=tmp[1,]
      data_tmp[p,]=tmp[2,]
    }
    data<-list()
    data$alphas<-alphas
    data$data<-data_tmp
    
    
    
  }
  if(distance_measure=="Bayes_estimated_priors"){
    
    alphas<-matrix(0, nrow=length(profile), ncol=length(subset) )
    priorgammas<-list()
    
    for(p in 1:length(profile)){
      print(p)
      alphas[p,]=apply( profile[[p]][,subset], 2, function (x) sum(x)/sum(summary[,p]) )
      #fit the gamma distribution
      
      zero_ind=which(alphas[p,]==0)
      
      if(length(zero_ind)!=0){ #some of alphas are zero
        
        priorgammas[[p]]=try(fitdistr(alphas[p,-zero_ind],"gamma")$estimate) #remove zero alphas
        if(class(priorgammas[[p]])=="try-error"){
          #print("err")
          priorgammas[[p]]=try(fitdistr(alphas[p, -zero_ind]/1000,"gamma")$estimate)
        }
      }
      else{ 
        priorgammas[[p]]=try(fitdistr(alphas[p,],"gamma")$estimate)
        if(class(priorgammas[[p]])=="try-error"){
          #print("err")
          priorgammas[[p]]=try(fitdistr(alphas[p,]/1000,"gamma")$estimate)
        }
        
        
      }
      
    }
    
    #compute the marginal likelihood
    marg_lh=matrix(0, nrow=length(profile), ncol=length(subset))
    for(p in 1:length(profile)){
      
      
      marg_lh[p,]=apply( profile[[p]][,subset], 2, FUN=marginal_likelihood, summary_p=summary[,p], A=priorgammas[[p]]["shape"], B=priorgammas[[p]]["rate"])
      
    }
    
    data=list()
    data$data=marg_lh
    data$priorgammas=priorgammas
    if(learn_alpha_prior==TRUE){
      data$alphas=alphas
      
    }
    
  }
  
  data
}


compute_distance_two_negatives<-function(profile, subset, summary_pos, summary_neg_promoters, summary_neg_random, fn=fn, distance_measure,
                                         learn_alpha_prior=TRUE, priorgammas_pos=NULL, priorgammas_neg_promoters=NULL, priorgammas_neg_random=NULL){
  
  
  if(distance_measure=="ML"){
    data_tmp<-matrix(0, nrow=3*length(profile), ncol=length(subset) ) #15 x N
    alpha_positives<-matrix(0, nrow=length(profile), ncol=length(subset) )  #15 x N
    alpha_negatives_promoters<-matrix(0, nrow=length(profile), ncol=length(subset))
    alpha_negatives_random<-matrix(0, nrow=length(profile), ncol=length(subset))
    for(p in 1:length(profile)){
      
      if(length(subset)!=1){
        tmp=apply( profile[[p]][,subset], 2, function (x) fn(x,summary_pos[,p]))
        alpha_positives[p,]=tmp[1,]
        data_tmp[p,]=tmp[2,]
      }
      else{
        tmp=fn(profile[[p]],summary_pos[,p])
        alpha_positives[p,]=tmp[1]
        data_tmp[p,]=tmp[2]
      }
      
    }
    
    for(p in 1:length(profile)){
      if(length(subset)!=1){
        tmp=apply( profile[[p]][,subset], 2, function (x) fn(x,summary_neg_promoters[,p]))
        alpha_negatives_promoters[p,]=tmp[1,]
        data_tmp[length(profile)+p,]=tmp[2,]
      }
      else{
        tmp=fn(profile[[p]],summary_neg_promoters[,p])
        alpha_negatives_promoters[p,]=tmp[1]
        data_tmp[length(profile)+p,]=tmp[2]
        
      }
      
    }
    
    for(p in 1:length(profile)){
      if(length(subset)!=1){
        tmp=apply( profile[[p]][,subset], 2, function (x) fn(x,summary_neg_random[,p]))
        alpha_negatives_random[p,]=tmp[1,]
        data_tmp[2*length(profile)+p,]=tmp[2,]
      }
      else{
        tmp=fn(profile[[p]],summary_neg_random[,p])
        alpha_negatives_random[p,]=tmp[1]
        data_tmp[2*length(profile)+p,]=tmp[2]
      }
      
    }
    
    
    
    data<-list()
    data$alpha_positives<-alpha_positives
    data$alpha_negatives_promoters<-alpha_negatives_promoters
    data$alpha_negatives_random<-alpha_negatives_random
    data$data<-data_tmp
    
    
    
  }
  if(distance_measure=="Bayes_estimated_priors"){
    
    if(learn_alpha_prior==TRUE){
      
      
      data_pos=learn_alpha_and_fit_gamma(profile=profile, subset=subset, summary=summary_pos)
      
      data_neg_promoters=learn_alpha_and_fit_gamma(profile=profile, subset=subset, summary=summary_neg_promoters)
      
      data_neg_random=learn_alpha_and_fit_gamma(profile=profile, subset=subset, summary=summary_neg_random)
      
      priorgammas_pos=data_pos$priorgammas    
      priorgammas_neg_promoters=data_neg_promoters$priorgammas
      priorgammas_neg_random=data_neg_random$priorgammas
      
      
    }
    
    #compute the marginal likelihood
    marg_lh=matrix(0, nrow=3*length(profile), ncol=length(subset))
    for(p in 1:length(profile)){
      
      if(length(subset)!=1){
        marg_lh[p,]=apply( profile[[p]][,subset], 2, FUN=marginal_likelihood, summary_p=summary_pos[,p], A=priorgammas_pos[[p]]["shape"], B=priorgammas_pos[[p]]["rate"])
      }
      else{
        marg_lh[p,]=marginal_likelihood(profile[[p]], summary_p=summary_pos[,p], A=priorgammas_pos[[p]]["shape"], B=priorgammas_pos[[p]]["rate"])
        
      }
    }
    
    for(p in 1:length(profile)){
      
      if(length(subset)!=1){
        marg_lh[length(profile)+p,]=apply( profile[[p]][,subset], 2, 
                                           FUN=marginal_likelihood, summary_p=summary_neg_promoters[,p], 
                                           A=priorgammas_neg_promoters[[p]]["shape"], B=priorgammas_neg_promoters[[p]]["rate"])
      }
      else{
        marg_lh[length(profile)+p,]=marginal_likelihood(profile[[p]], summary_p=summary_neg_promoters[,p], 
                                                        A=priorgammas_neg_promoters[[p]]["shape"], B=priorgammas_neg_promoters[[p]]["rate"])
        
      }
    }
    
    
    for(p in 1:length(profile)){
      if(length(subset)!=1){
        marg_lh[2*length(profile)+p,]=apply( profile[[p]][,subset], 2, 
                                             FUN=marginal_likelihood, summary_p=summary_neg_random[,p], 
                                             A=priorgammas_neg_random[[p]]["shape"], B=priorgammas_neg_random[[p]]["rate"])
      }
      else{
        marg_lh[2*length(profile)+p,]=marginal_likelihood(profile[[p]], summary_p=summary_neg_random[,p], 
                                                          A=priorgammas_neg_random[[p]]["shape"], B=priorgammas_neg_random[[p]]["rate"])
        
      }
    }
    
    
    data=list()
    data$data=marg_lh
    data$priorgammas_pos=priorgammas_pos
    data$priorgammas_neg_promoters=priorgammas_neg_promoters
    data$priorgammas_neg_random=priorgammas_neg_random
    if(learn_alpha_prior==TRUE){
      data$alpha_positives<-data_pos$alpha
      data$alpha_negatives_promoters<-data_neg_promoters$alpha
      data$alpha_negatives_random<-data_neg_random$alpha
    }
    
  }
  
  
  data
}


learn_alpha_and_fit_gamma<-function(profile, subset, summary){
  
  alpha<-matrix(0, nrow=length(profile), ncol=length(subset) )
  priorgammas<-list()
  
  for(p in 1:length(profile)){
    
    alpha[p,]=apply( profile[[p]][,subset], 2, function (x) sum(x)/sum(summary[,p]) )
    #fit the gamma distribution
    
    zero_ind=which(alpha[p,]==0)
    
    if(length(zero_ind)!=0){ #some of alphas are zero
      
      priorgammas[[p]]=try(fitdistr(alpha[p,-zero_ind],"gamma")$estimate) #remove zero alphas
      if(class(priorgammas[[p]])=="try-error"){
        #print("err")
        priorgammas[[p]]=try(fitdistr(alpha[p, -zero_ind]/1000,"gamma")$estimate)
      }
    }
    else{ 
      priorgammas[[p]]=try(fitdistr(alpha[p,],"gamma")$estimate)
      if(class(priorgammas[[p]])=="try-error"){
        #print("err")
        priorgammas[[p]]=try(fitdistr(alpha[p,]/1000,"gamma")$estimate)
      }
      
      
    }
    
  } 
  
  data<-list()
  data$alpha=alpha
  data$priorgammas=priorgammas
  data
  
}

marginal_likelihood<-function(profile_p, summary_p, A, B){
  value=A*log(B)+sum( profile_p * log(summary_p) )-( A + sum(profile_p) ) * log( B + sum(summary_p) )+lgamma( A + sum (profile_p) )-lgamma( A )-sum ( lgamma( profile_p + 1) )
  
  value
}

compute_Gamma_posterior<-function(profiles_p, summary_p, a, b){
  
  A=a+sum( profiles_p  )
  B=b+ncol(profiles_p)*sum(summary_p)
  
  c(A,B)
  
}




compute_ML<-function(data_profile,summary){
  #data_profile is a vector
  #summary is a vector
  alpha=sum(data_profile)/sum(summary)
  
  lh=sum( dpois(x=data_profile ,lambda=alpha*summary, log=TRUE) )
  
  rt=c(alpha, lh)
  names(rt)=c("alpha","lh")
  rt
}





compute_Bayes_estimated_prior<-function(){
  #learn the prior from the data, compute the marginal likelihood for all samples
  
}


compute_Bayes_fixed_prior<-function(){
  
}





createCrossValidationGroups<- function(N_pos, N_promoter, N_random, k){
  
  N_neg=N_promoter+N_random
  pos_ind_rand=sample(seq(1,N_pos,1),N_pos)
  
  neg_ind_rand=sample(seq(1,N_neg,1),N_neg)
  
  pos_ind_cv<-list()
  neg_ind_cv<-list() 
  
  #divide the data into k parts
  for(i in 1:k){
    
    pos_ind_cv[[i]]<-pos_ind_rand[((i-1)*(length(pos_ind_rand)/k)+1):(i*( length(pos_ind_rand) /k))]
    neg_ind_cv[[i]]<-neg_ind_rand[((i-1)*( length(neg_ind_rand) /k)+1):(i*( length(neg_ind_rand) /k))]
    
    
    
  }    
  
  pos_groups_cv<-list()
  neg_groups_cv<-list()
  
  #cv groups
  
  for(i in 1:k){
    pos_test_ind<-pos_ind_cv[[i]] #the test data
    pos_train_ind<-c(pos_ind_cv[-i], recursive=TRUE) #the rest is train data
    
    neg_test_ind<-neg_ind_cv[[i]]
    neg_train_ind<-c(neg_ind_cv[-i],recursive=TRUE)
    
    tmp<-list()
    tmp$test<-pos_test_ind
    tmp$train<-pos_train_ind
    pos_groups_cv[[i]]<-tmp
    
    rm(tmp)
    tmp<-list()
    tmp$test<-neg_test_ind
    tmp$train<-neg_train_ind
    neg_groups_cv[[i]]<-tmp
    
  }
  
  neg_groups_cv_promoters<-list()
  neg_groups_cv_random<-list()
  
  for(i in 1:length(neg_groups_cv)){
    rm(tmp)
    tmp<-list()
    tmp$train=neg_groups_cv[[i]]$train[which(neg_groups_cv[[i]]$train<=N_promoter)]
    tmp$test=neg_groups_cv[[i]]$test[which(neg_groups_cv[[i]]$test<=N_promoter)]
    neg_groups_cv_promoters[[i]]=tmp
  }
  
  for(i in 1:length(neg_groups_cv)){
    rm(tmp)
    tmp<-list()
    tmp$train=neg_groups_cv[[i]]$train[which(neg_groups_cv[[i]]$train>N_promoter)]
    tmp$test=neg_groups_cv[[i]]$test[which(neg_groups_cv[[i]]$test>N_promoter)]
    neg_groups_cv_random[[i]]=tmp
  }
  
  
  
  returnList<-list()
  returnList$pos=pos_groups_cv
  returnList$neg_promoter=neg_groups_cv_promoters
  returnList$neg_random=neg_groups_cv_random
  returnList
}