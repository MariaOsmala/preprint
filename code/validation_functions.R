enhancer_grouping_choose_type<-function(predictions_full, filename, type, 
                                        threshold, enhancer_separation, overlap){
  #type is maxscore or middle or multi
  #enhancer_separation is only used when type="multi", the minimal distance between adjacent enhancers
  enhancer_ind=which(predictions_full$enhancer_score>threshold) #548150
  window=unique(width(ranges(predictions_full)))
  #the distance between subsequent enhancers in predictions_full
  enhancers_subseq_index=enhancer_ind[2:length(enhancer_ind)]-enhancer_ind[1:(length(enhancer_ind)-1)]
  if(enhancer_ind[2]==enhancer_ind[1]+1){   #if the first two values are the same, both ones -> subsequent
    enhancers_subseq_index=c(0,enhancers_subseq_index) 
  }
  
  
  test=rle2_better_works(enhancers_subseq_index,indices = TRUE, nmax_input=10000)
  
  i=2
  while(test[nrow(test),3]!=length(enhancers_subseq_index)){
    test=rle2_better_works(enhancers_subseq_index,indices = TRUE, nmax_input=i*10000)
    i=i+1
  } 
  
  
  #each row must have 1 either in the first or the fourth column
  problem=which(test[,1]!=1 & test[,4]!=1) 
  
  while(length(problem)!=0 ){
    
    test1=test[1:(problem[1]-1),] #from the beginning to the problem
    
    #divide the problematic regions to individual regions
    test3=as.matrix(data.frame(values=test[problem[1],1]:(test[problem[1],1]+test[problem[1],4]-1), 
                               starts=test[problem[1],2]:test[problem[1],3], 
                               ends=test[problem[1],2]:test[problem[1],3], lengths=rep(1,test[problem[1],4]) ))
    
    if(problem[1]!=nrow(test)){
      test2=test[(problem[1]+1):nrow(test),] #from the problem to the end
      
      test=rbind(test1, test3, test2)
    }else{
      test=rbind(test1, test3)
    }
    problem=which(test[,1]!=1 & test[,4]!=1)
    
    
  }
  
  
  
  #these are the regions with subsequent enhancer predictions
  oneind=which(test[,1]==1) #what are these in the original indexing
  
  if(enhancer_ind[2]== (enhancer_ind[1]+1)){ #immediately in the beginnin these are subsequent enhancers
    grouped_start=test[oneind-1,2]
    grouped_end=test[oneind,3]
  }else{				#there is an individual region in the beginning
    grouped_start=test[oneind-1,2]+1
    grouped_end=test[oneind,3]+1
  }
  
  #these are the indexes in enhancer_ind for wide subsequent regions
  grouped<-data.frame(start=grouped_start, end=grouped_end)
  
  #remove grouped from the enhancer list
  seq_vector <- function(x, c1) {
    #print(x)
    seq(x[1], x[2], c1)
  }
  byy=1
  
  #indeces of enhancer_ind in a list
  tmp=apply(as.matrix(grouped), 1, seq_vector,  c1 = 1) #subsequent enhancer regions, indexes in enhancers_ind
  
  #same but unlisted	
  tmp_all=unlist(tmp)
  
  #indiv_enhancers=enhancer_ind[-tmp_all] #617
  
  
  if(type=="middle"){
    
    grouped_middle<-GRangesList()
    max_scores<-vector() 
    for(i in 1:length(tmp)){
      #print(i)
      
      grouped_middle[[i]]=resize(reduce(predictions_full[enhancer_ind[tmp[[i]] ] ]), window, fix="center")
      max_scores[i]=max(predictions_full$enhancer_score[ enhancer_ind[ tmp[[i]]] ])
      
    }
    
    predictions_middle_full=predictions_full[-enhancer_ind[tmp_all]]
    
    test=unlist(grouped_middle) #some wide subsequent enhancer regions are missing, like
    elementMetadata(test)=data.frame(label=rep("1", length(test)), enhancer_score=max_scores)
    predictions_middle_full=c(predictions_middle_full, test)
    result=predictions_middle_full
    
    
  }
  else if(type=="maxscore"){ #experimental
    maxscore_locations<-vector()
    smallscore_locations<-list()
    max_scores<-vector()
    max_scores_length<-vector()
    for(i in 1:length(tmp)){
      #print(i)
      max_scores[i]=max(predictions_full$enhancer_score[ enhancer_ind[ tmp[[i]]] ])
      max_scores_length[i]= length( which(predictions_full$enhancer_score[ enhancer_ind[ tmp[[i]] ] ]==max_scores[i]) )
      #the max score can be in multiple places, sample randomly from these
      if(max_scores_length[i]!=1){
        maxscore_locations[i]=sample( enhancer_ind[tmp[[i]]][which(predictions_full$enhancer_score[enhancer_ind[ tmp[[i]] ] ]==max_scores[i])] ,1)
      }else{
        maxscore_locations[i]=enhancer_ind[tmp[[i]]][which(predictions_full$enhancer_score[enhancer_ind[ tmp[[i]] ] ]==max_scores[i])] 
      }
      smallscore_locations[[i]]=enhancer_ind[tmp[[i]] ][-which( enhancer_ind[tmp[[i]] ]==  maxscore_locations[i])]
      
      
      
    }
    predictions_maxscore_full=predictions_full[-unlist(smallscore_locations)]
    
    result=predictions_maxscore_full
    
  }
  else if(type=="multi"){ #experimental
    
    rm(multi_scores)
    rm(multi_locations)
    multi_scores<-vector()
    multi_locations_start<-vector()
    multi_locations_end<-vector()
    
    peak_dist=enhancer_separation/overlap #minimal distance between individual enhancers OVERLAP or BIN_SIZE
    #prob_dist #the prediction scores of a single enhancer region composed of adjacent regions having an enhancer score larger than threshold
    for(tmp_ind in 1:length(tmp)){
      
      
      prob_dist=predictions_full[enhancer_ind[tmp[[tmp_ind]]]]$enhancer_score
      
      region_start=as.data.frame(ranges(predictions_full[enhancer_ind[tmp[[tmp_ind]]]]))$start[1]
      
      rm(indiv_enhancers_pks)
      rm(indiv_enhancers_locs)  
      indiv_enhancers_pks<-vector()
      indiv_enhancers_locs<-vector()
      
      
      #prob_dist=predictions_full[enhancer_ind[tmp[[211]]]]$enhancer_score #vector of length 147, this many subsequent enhancer predictions
      
      if(length(prob_dist)<=peak_dist){
        #the length of a wide enhancer region is equal or less than 2000 bp
        #finds the maximum value -> pks
        #locs is the first occurrence of the maximal value
        pks=max(prob_dist);
        locs=which(prob_dist==pks)[1]
        indiv_enhancers_pks=c(indiv_enhancers_pks, pks)
        indiv_enhancers_locs=c(indiv_enhancers_locs, locs)  
        
      }
      else{
        #the length of a wider enhancer region is more than 2000 bp	
        
        pks=max(prob_dist);
        locs=which(prob_dist==pks)[1] #the first occurrence of the maximal value
        indiv_enhancers_pks=c(indiv_enhancers_pks, pks)
        indiv_enhancers_locs=c(indiv_enhancers_locs, locs)  
        first_vector_end=locs-peak_dist
        second_vector_start=locs+peak_dist
        rm(vectors)
        vectors<-list()
        if(first_vector_end>0){
          rm(vector)
          vector<-list()
          vector[[1]] <-prob_dist[1:first_vector_end]
          vector[[2]]<-1 # index of start in the original vector
          vectors[[length(vectors)+1]]<-vector
        }
        if(second_vector_start<= length(prob_dist)){
          vector<-list()
          vector[[1]] <-prob_dist[second_vector_start:length(prob_dist)]
          vector[[2]]<-second_vector_start	
          
          vectors[[length(vectors)+1]]<-vector
        }
        
        compute_vector_length<-function(tmp){
          length(tmp[[1]])
        }
        
        #iiii=1	
        while(length(vectors)!=0){
          #print(iiii)
          #iiii=iiii+1
          #are there vectors with length less or equal to peak_dist
          
          #which of the vectors is shorter or equal to 2000 bp
          peak_dist_length=which(unlist(lapply(vectors, compute_vector_length)) <= peak_dist)
          if(length(peak_dist_length)>0){
            for(i in peak_dist_length){
              pks=max(vectors[[i]][[1]])
              locs=which(vectors[[i]][[1]]==pks )[1]
              indiv_enhancers_pks=c(indiv_enhancers_pks, pks)
              indiv_enhancers_locs=c(indiv_enhancers_locs, vectors[[i]][[2]]+locs-1)
              
              
            }
            #remove peak_dist_length corresponding vectors, shorter or equal to 2000 bp
            if(length(peak_dist_length)<length(vectors) ){
              vectors<-vectors[-peak_dist_length]
            }
            else{
              vectors<-NULL
            }		
            
          }
          
          
          if(length(vectors)!=0){
            #which of the vectors has length larger than 2000 bp	
            peak_dist_length=which(unlist(lapply(vectors, compute_vector_length)) > peak_dist)
            
          }
          else{
            peak_dist_length=NULL
          }
          
          if(length(peak_dist_length)!=0){
            
            for(i in peak_dist_length){
              pks=max(vectors[[i]][[1]])
              locs=which(vectors[[i]][[1]]==pks )[1]
              indiv_enhancers_pks=c(indiv_enhancers_pks, pks)
              indiv_enhancers_locs=c(indiv_enhancers_locs, vectors[[i]][[2]]+locs-1)
              first_vector_end=locs-peak_dist
              second_vector_start=locs+peak_dist
              
              if(first_vector_end>0){
                rm(vector)
                vector<-list()
                vector[[1]] <-vectors[[i]][[1]][1:first_vector_end]
                vector[[2]]<-vectors[[i]][[2]]
                vectors[[length(vectors)+1]]<-vector
              }
              if(second_vector_start<= length(vectors[[i]][[1]]) ){
                rm(vector)	        	
                vector<-list()
                vector[[1]] <-vectors[[i]][[1]][second_vector_start:length(vectors[[i]][[1]])]
                
                vector[[2]]<-vectors[[i]][[2]]+second_vector_start-1	
                
                vectors[[length(vectors)+1]]<-vector
              }
              
              
            } #close for
            
            #remove those vectors that were handled, longer than 2000 bp
            if(length(peak_dist_length)<length(vectors) ){ #1 99 1 78 99 164 	# 59 164 38 185 26 197
              vectors<-vectors[-peak_dist_length]
            }
            else{
              vectors<-NULL
            }	
            
            
            
          } #close if
          
          
          
          
          
          
        } #while closes here
        
        
        
      } #else closes here
      
      multi_scores<-c(multi_scores, indiv_enhancers_pks)
      multi_locations_start<-c(multi_locations_start, region_start+(indiv_enhancers_locs-1)*overlap)  #overlap, who much the region is shifted
      multi_locations_end<-c(multi_locations_end, region_start+(indiv_enhancers_locs-1)*overlap+window-1)
      
      
    } #tmp closes here	
    
    
    predictions_multi_full=predictions_full[-enhancer_ind[tmp_all]]
    
    
    test=GRanges( seqnames = Rle( rep(chrlist, length(multi_locations_start)) , rep(1, length(multi_locations_start) ) ), 
                  ranges = IRanges(start=multi_locations_start, end=multi_locations_end), 
                  strand = Rle( strand( rep( '+',length(multi_locations_start) ) ), rep(1,length(multi_locations_start) ) ), 
                  label=rep(1,length(multi_locations_start)), enhancer_score=multi_scores)
    
    
    
    predictions_multi_full=c(predictions_multi_full, test)
    result=predictions_multi_full
    
    
  }	#multi closes here
  else{
    print("error")
  }
  
  
  
  result
}



#x=enhancers_subseq_index
#indices=TRUE
#return.list=FALSE
rle2_better<-function (x, indices = FALSE, return.list = FALSE, nmax_input=10000) 
{
  if (is.matrix(x) | is.data.frame(x)) { #FALSE
    if (ncol(x) > 1) {
      stop("For x= option, please enter vector or single-column matrix or data frame (can be numeric or character)")
    }
  }
  if (!is.logical(indices)) { #FALSE
    stop("For indices= option, please enter TRUE or FALSE")
  }
  if (!is.logical(return.list)) { #FALSE
    stop("For return.list= option, please enter TRUE or FALSE")
  }
  if (is.numeric(x)) { #TRUE
    class_x = 1
  }
  else if (is.character(x)) { #FALSE
    class_x = 2
  }
  else {
    stop("For x= option, please enter vector or single-column matrix or data frame (can be numeric or character)")
  }
  indices = ifelse(indices, 1, 0) #1
  length_x = length(x) #46759
  if (length_x <= 10000) { #FALSE
    if (class_x == 1) {
      out = rle2.num(x = x, n = length_x, nmax = -1, indices = indices)
    }
    else {
      out = rle2.char(x = x, n = length_x, nmax = -1, indices = indices)
    }
  }
  else {
    end_partial = ceiling(length_x/1000) #47
    if (class_x == 1) {
      partial = rle2.num(x = x[1:end_partial], n = end_partial, 
                         nmax = -1, indices = 0)
    }
    else {
      partial = rle2.char(x = x[1:end_partial], n = end_partial, 
                          nmax = -1, indices = 0)
    }
    rows = nrow(partial) #6
    if (rows/end_partial > 0.5) { #FALSE
      if (class_x == 1) {
        out = rle2.num(x = x, n = length_x, nmax = -1, 
                       indices = indices)
      }
      else {
        out = rle2.char(x = x, n = length_x, nmax = -1, 
                        indices = indices)
      }
    }
    else {
      nmax = nmax_input
      
      repeat {
        if (class_x == 1) {
          out = rle2.num(x = x, n = length_x, nmax = nmax, 
                         indices = indices)
        }
        else {
          out = rle2.char(x = x, n = length_x, nmax = nmax, 
                          indices = indices)
        }
        if (nrow(out) < nmax) {
          break
        }
        else {
          nmax = nmax * 10
        }
      }
    }
  }
  if (class_x == 2) {
    if (indices == 0) {
      out = cbind(as.data.frame(x = cbind(out[, 1]), stringsAsFactors = FALSE), 
                  as.numeric(out[, 2]))
    }
    else {
      out = cbind(as.data.frame(x = cbind(out[, 1]), stringsAsFactors = FALSE), 
                  as.numeric(out[, 2]), as.numeric(out[, 3]), as.numeric(out[, 
                                                                             4]))
    }
  }
  if (indices == 0) {
    colnames(out) = c("values", "lengths")
  }
  else {
    colnames(out) = c("values", "starts", "stops", "lengths")
  }
  if (return.list == TRUE) {
    if (indices == 0) {
      out = list(values = out[, 1], lengths = out[, 2])
    }
    else {
      out = list(values = out[, 1], starts = out[, 2], 
                 stops = out[, 3], lengths = out[, 4])
    }
  }
  return(out)
}

rle2_better_works<-function (x, indices = FALSE, return.list = FALSE, nmax_input=10000) 
{
  if (is.matrix(x) | is.data.frame(x)) { #FALSE
    if (ncol(x) > 1) {
      stop("For x= option, please enter vector or single-column matrix or data frame (can be numeric or character)")
    }
  }
  if (!is.logical(indices)) { #FALSE
    stop("For indices= option, please enter TRUE or FALSE")
  }
  if (!is.logical(return.list)) { #FALSE
    stop("For return.list= option, please enter TRUE or FALSE")
  }
  if (is.numeric(x)) { #TRUE
    class_x = 1
  }
  else if (is.character(x)) { #FALSE
    class_x = 2
  }
  else {
    stop("For x= option, please enter vector or single-column matrix or data frame (can be numeric or character)")
  }
  indices = ifelse(indices, 1, 0) #1
  length_x = length(x) #46759
  if (length_x <= 10000) { #FALSE
    if (class_x == 1) {
      out = .Call("accelerometry_rle2_num",PACKAGE="accelerometry",x = x, n = length_x, nmax = -1, indices = indices)
    }
    else {
      out = .Call("accelerometry_rle2_char",PACKAGE="accelerometry",x = x, n = length_x, nmax = -1, indices = indices)
      
    }
  }
  else {
    end_partial = ceiling(length_x/1000) #47
    if (class_x == 1) {
      partial = .Call("accelerometry_rle2_num",PACKAGE="accelerometry",x = x[1:end_partial], n = end_partial, nmax = -1, indices = 0)
      
      
    }
    else {
      partial = .Call("accelerometry_rle2_char",PACKAGE="accelerometry",x = x[1:end_partial], n = end_partial, nmax = -1, indices = 0)
      
    }
    rows = nrow(partial) #6
    if (rows/end_partial > 0.5) { #FALSE
      if (class_x == 1) {
        out = .Call("accelerometry_rle2_num",PACKAGE="accelerometry",x = x, n = length_x, nmax = -1, indices = indices) 
        
        
      }
      else {
        out = .Call("accelerometry_rle2_char",PACKAGE="accelerometry",x = x, n = length_x, nmax = -1, indices = indices)
        
      }
    }
    else {
      nmax = nmax_input
      
      repeat {
        if (class_x == 1) {
          out = .Call("accelerometry_rle2_num",PACKAGE="accelerometry",x = x, n = length_x, nmax = nmax, indices = indices) 
          
        }
        else {
          out = .Call("accelerometry_rle2_char",PACKAGE="accelerometry",x = x, n = length_x, nmax = nmax, indices = indices)
          
        }
        if (nrow(out) < nmax) {
          break
        }
        else {
          nmax = nmax * 10
        }
      }
    }
  }
  if (class_x == 2) {
    if (indices == 0) {
      out = cbind(as.data.frame(x = cbind(out[, 1]), stringsAsFactors = FALSE), 
                  as.numeric(out[, 2]))
    }
    else {
      out = cbind(as.data.frame(x = cbind(out[, 1]), stringsAsFactors = FALSE), 
                  as.numeric(out[, 2]), as.numeric(out[, 3]), as.numeric(out[, 
                                                                             4]))
    }
  }
  if (indices == 0) {
    colnames(out) = c("values", "lengths")
  }
  else {
    colnames(out) = c("values", "starts", "stops", "lengths")
  }
  if (return.list == TRUE) {
    if (indices == 0) {
      out = list(values = out[, 1], lengths = out[, 2])
    }
    else {
      out = list(values = out[, 1], starts = out[, 2], 
                 stops = out[, 3], lengths = out[, 4])
    }
  }
  return(out)
}

findOverlaps_select_nearest <- function(query, subject, minoverlap=1, type="any"){
  
  mymethod2RFECS<-findOverlaps(query, subject, minoverlap=minoverlap, type="any")
  print(length(mymethod2RFECS)  )
  
  rle_from=rle(from(mymethod2RFECS))
  multiple=rle_from$values[which(rle_from$lengths > 1)]
  rm(remove_ind_vector)
  remove_ind_vector<-c()
  max_ind_vector<-c()
  
  for(i in 1:length(multiple)){
    my_enhancer=query[multiple[i]]
    RFECS_enhancers_set=subject[to(mymethod2RFECS)[which( from(mymethod2RFECS) == multiple[i])]]
    mymethod2RFECS_overlap<-c()
    for(j in 1:length(RFECS_enhancers_set)){
      mymethod2RFECS_overlap[j] <- width(pintersect(my_enhancer, RFECS_enhancers_set[j]))
    }
    max_ind=which.max(mymethod2RFECS_overlap)
    #max_ind_vector[i]=max_ind
    #print(paste("max_ind: ",max_ind,sep=""))
    #remove all other overlaps except the maximum
    
    remove_ind=1:length(mymethod2RFECS_overlap)
    remove_ind=remove_ind[-max_ind]
    
    remove_ind_vector=c(remove_ind_vector, which( from(mymethod2RFECS) == multiple[i])[remove_ind] )
    
    
    
  }
  
  print(length(remove_ind_vector))
  mymethod2RFECS<- mymethod2RFECS[-remove_ind_vector]
  
  #Another way around
  rle_to=rle(to(mymethod2RFECS))
  multiple=rle_to$values[which(rle_to$lengths > 1)]
  rm(remove_ind_vector)
  remove_ind_vector<-c()
  max_ind_vector<-c()
  
  for(i in 1:length(multiple)){
    RFECS_enhancer=subject[multiple[i]]
    mymethod_enhancers_set=query[from(mymethod2RFECS)[which( to(mymethod2RFECS) == multiple[i])]]
    mymethod2RFECS_overlap<-c()
    for(j in 1:length(RFECS_enhancers_set)){
      mymethod2RFECS_overlap[j] <- width(pintersect(RFECS_enhancer, mymethod_enhancers_set[j]))
    }
    max_ind=which.max(mymethod2RFECS_overlap)
    #max_ind_vector[i]=max_ind
    #print(paste("max_ind: ",max_ind,sep=""))
    #remove all other overlaps except the maximum
    
    remove_ind=1:length(mymethod2RFECS_overlap)
    remove_ind=remove_ind[-max_ind]
    
    remove_ind_vector=c(remove_ind_vector, which( to(mymethod2RFECS) == multiple[i])[remove_ind] )
    
    
    
  }
  
  
  print(length(remove_ind_vector))
  mymethod2RFECS<- mymethod2RFECS[-remove_ind_vector]
  
  
  
  return(mymethod2RFECS)
  
}