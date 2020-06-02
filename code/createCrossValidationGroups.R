createCrossValidationGroups<- function(N_pos, N_promoter, N_random, k){
  
  #there should be equal amount of promoters and random regions in each CV-set
  #N_neg=N_promoter+N_random
  
  #First positive
  pos_ind_rand=sample(seq(1,N_pos,1),N_pos) #1000 in random order
  #neg_ind_rand=sample(seq(1,N_neg,1),N_neg) #2000/3000 in random order
  
  pos_ind_cv<-list()
  #neg_ind_cv<-list() 
  
  #divide the data into k parts, first the positive
  for(i in 1:k){
    pos_ind_cv[[i]]<-pos_ind_rand[((i-1)*(length(pos_ind_rand)/k)+1):(i*( length(pos_ind_rand) /k))] # length 200, number of pos/k=1000/5
    #neg_ind_cv[[i]]<-neg_ind_rand[((i-1)*( length(neg_ind_rand) /k)+1):(i*( length(neg_ind_rand) /k))] #length 600, number of neg/k=3000/5
  }    
  pos_groups_cv<-list()
  for(i in 1:k){
    pos_test_ind<-pos_ind_cv[[i]] #the test data 200
    pos_train_ind<-c(pos_ind_cv[-i], recursive=TRUE) #the rest is train data, 800
    
    #neg_test_ind<-neg_ind_cv[[i]] #600
    #neg_train_ind<-c(neg_ind_cv[-i],recursive=TRUE) #2400
    
    tmp<-list()
    tmp$test<-pos_test_ind #200
    tmp$train<-pos_train_ind #800
    pos_groups_cv[[i]]<-tmp #list of 2
    
    #rm(tmp)
    #tmp<-list()
    #tmp$test<-neg_test_ind
    #tmp$train<-neg_train_ind
    #neg_groups_cv[[i]]<-tmp
    
  } 
  
  ######################Promoters####################################
  
  promoter_ind_rand=sample(seq(1,N_promoter,1),N_promoter) #1000 in random order
  promoter_ind_cv<-list()
  
  #divide the promoters into k parts
  for(i in 1:k){
    promoter_ind_cv[[i]]<-promoter_ind_rand[((i-1)*(length(promoter_ind_rand)/k)+1):(i*( length(promoter_ind_rand) /k))] # length 200, number of pos/k=1000/5
    
  }    
  
  promoter_groups_cv<-list()
  for(i in 1:k){
    promoter_test_ind<-promoter_ind_cv[[i]] #the test data 200
    promoter_train_ind<-c(promoter_ind_cv[-i], recursive=TRUE) #the rest is train data, 800
    
    tmp<-list()
    tmp$test<-promoter_test_ind #200
    tmp$train<-promoter_train_ind #800
    promoter_groups_cv[[i]]<-tmp #list of 2
  } 
  
  
  ######################Random####################################
  
  random_ind_rand=sample(seq(1,N_random,1),N_random) #1000 or 2000 indexes in random order
  random_ind_cv<-list()
  
  #divide the promoters into k parts
  for(i in 1:k){
    random_ind_cv[[i]]<-random_ind_rand[((i-1)*(length(random_ind_rand)/k)+1):(i*( length(random_ind_rand) /k))] # length 200, number of pos/k=1000/5
    
  }    
  
  random_groups_cv<-list()
  for(i in 1:k){
    random_test_ind<-random_ind_cv[[i]] #the test data 200
    random_train_ind<-c(random_ind_cv[-i], recursive=TRUE) #the rest is train data, 800
    
    tmp<-list()
    tmp$test<-random_test_ind #200
    tmp$train<-random_train_ind #800
    random_groups_cv[[i]]<-tmp #list of 2
  } 
  
  
  
  
  returnList<-list()
  returnList$pos=pos_groups_cv
  returnList$neg_promoter=promoter_groups_cv
  returnList$neg_random=random_groups_cv
  returnList
}