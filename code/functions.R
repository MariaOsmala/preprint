create_promoter_list<-function(TSS_annotation, TSS_annotation_positive, DNase_peaks_file, N,
                               between_TSS_distance, N_initial_selection, 
                               remove_blacklist_regions=TRUE, ENCODE_blacklist){
  
  
  #N is the number of regions to be selected
  
  #take the direction into account or not
  
  #consider only those TSS that do not have any other TSS near-by
  
  #choose the TSS that is closests to the peak summit of DNase-peak
    
  
  GR_DNase=narrowPeak2GRanges(narrowPeakFile = DNase_peaks_file, asPeak=FALSE, elementMetadata_significance = "pValue")
  GR_DNase_peak=narrowPeak2GRanges(narrowPeakFile = DNase_peaks_file, asPeak=TRUE, elementMetadata_significance = "pValue")
  
  
  test=as.data.frame( distanceToNearest(GR_DNase_peak,TSS_annotation_positive))
  
  #sort by the distance
  test=test[order(test$distance, decreasing=FALSE),]
  
  test_subject=test$subjectHits[1:N_initial_selection]
  test_query=test$queryHits[1:N_initial_selection]
  
  TSS_with_DNase=TSS_annotation_positive[test_subject]
  DNase_signalValue=elementMetadata(GR_DNase[test_query])$peaks.signalValue
  DNase_pvalue=elementMetadata(GR_DNase[test_query])$peaks.pValue
  
  elementMetadata(TSS_with_DNase)=data.frame(signalValue=DNase_signalValue, pValue=DNase_pvalue)
  
  
  
  TSS_strand=strand(TSS_annotation[test_subject])
  
  test=as.matrix(findOverlaps(TSS_with_DNase, GR_DNase))
  
  TSS_with_DNase=TSS_with_DNase[unique(test[,1])]
  
  TSS_strand=TSS_strand[unique(test[,1])]
  
  #test=as.matrix(findOverlaps(GR_Gencode_all_positive, GR_DNase))
  
  #TSS_with_DNase=GR_Gencode_all_positive[test[,1]]
  
  #TSS_strand=strand(GR_Gencode_all[test[,2]])
  
  
  
  i=1
  while(i<=length(TSS_with_DNase)){
    test2=distance(TSS_with_DNase[i], TSS_with_DNase) 
    too_close=which(test2<=between_TSS_distance & test2!=0)
    
    if(length(too_close)!=0){
      TSS_with_DNase=TSS_with_DNase[-too_close]
      TSS_strand=TSS_strand[-too_close]
      
    }
    
    i=i+1
    
  }
  
  
  
  
  #elementMetadata(TSS_with_DNase)=elementMetadata(GR_DNase[test[,2]])[,paste("elementMetadata.",order,sep="")]
  #order based on the "order"
  
  order_ind=order(TSS_with_DNase$pValue, decreasing=TRUE)
  TSS_with_DNase=TSS_with_DNase[order_ind]
  TSS_strand=TSS_strand[order_ind]
  order_ind=order(TSS_with_DNase$signalValue, decreasing=TRUE)
  TSS_with_DNase=TSS_with_DNase[order(TSS_with_DNase$signalValue, decreasing=TRUE)]
  TSS_strand=TSS_strand[order_ind]
    
 
  
  regions_wider=resize(TSS_with_DNase, width=2000, fix="center" )
  #what about the blacklist regions, are these zero or 1-based
  #blacklist regions are 1-based
  if(remove_blacklist_regions==TRUE){
    test=as.matrix(findOverlaps(TSS_with_DNase, ENCODE_blacklist))
    if(nrow(test)!=0){
      TSS_with_DNase=TSS_with_DNase[-unique(test[,1])]
      TSS_strand=TSS_strand[-unique(test[,1])]
    }
    
  }
  
  
  return_list<-list()
  return_list$promoters<-TSS_with_DNase
  return_list$strand=TSS_strand
  return_list
  
  
}




narrowPeak2GRanges<-function(narrowPeakFile, asPeak=FALSE, elementMetadata_significance="qValue"){
  
  
  #narrowPeak format == BED6+4 format, zero based
  library(BSgenome.Hsapiens.UCSC.hg19)
  human.chromlens = seqlengths(Hsapiens)
  
  peaks=read.table(narrowPeakFile, stringsAsFactors=FALSE)    
  names(peaks)<-c("chrom","chromStart", "chromEnd", "name", "score","strand", "signalValue", "pValue","qValue", "peak")
  
  if(asPeak==TRUE){
    #peak - Point-source called for this peak; 0-based offset from chromStart. Use -1 if no point-source called.
    peaksGRanges=GRanges(seqnames = Rle(peaks$chrom, rep(1, length(peaks$chrom)) ), 
                                ranges = IRanges(start=peaks$chromStart+peaks$peak+1, end=peaks$chromStart+peaks$peak+1),  
                                strand = Rle( strand( rep( '+',nrow(peaks) ) )) )
  }else{
    
    
    peaksGRanges=GRanges(seqnames = Rle(peaks$chrom, rep(1, length(peaks$chrom)) ), 
                         ranges = IRanges(start=peaks$chromStart+1, end=peaks$chromEnd),  
                         strand = Rle( strand( rep( '+',nrow(peaks) ) )) )
    
  }
  
  #for DNase-seq this is pValue, for p300 signalValue
  if(elementMetadata_significance=="qValue"){
      elementMetadata(peaksGRanges) <- data.frame(peaks$qValue, peaks$signalValue)
  }else{
    elementMetadata(peaksGRanges) <- data.frame(peaks$pValue, peaks$signalValue)
  }
  
  seqlengths(peaksGRanges)<-human.chromlens[seqnames(seqinfo(peaksGRanges))]

  peaksGRanges  
  
  
    
  
  
  
  
  
  
}



removeChrFromGRanges <- function(peaksGRanges, chr){
  if(length(which(seqnames(peaksGRanges)==chr))!=0){
    peaksGRanges=peaksGRanges[-which(seqnames(peaksGRanges)==chr)]
    
  }
  peaksGRanges
  
}



create_enhancer_list<-function(p300, DNase_peaks_file,  TSS_annotation, distance_to_promoters, remove_blacklist_regions=TRUE, ENCODE_blacklist){
    
   
    
    
    p300peaks_GRanges=narrowPeak2GRanges(p300, TRUE, "qValue")
   
    
  
   
    p300peaks_GRanges=removeChrFromGRanges(p300peaks_GRanges, "chrM")
  
  
    
    close_to_TSS=which(abs(distToTss(p300peaks_GRanges,TSS_annotation) )< distance_to_promoters) 
    GR_Enhancers=p300peaks_GRanges[-close_to_TSS]
    
   
   
    
    DNase_peaks_GRanges=narrowPeak2GRanges(DNase_peaks_file, FALSE, 'pValue')
    
   
   
    
    
    test=as.matrix(findOverlaps(GR_Enhancers, DNase_peaks_GRanges)) 
        

    GR_Enhancers=GR_Enhancers[test[,1]]
    GR_Enhancers=GR_Enhancers[order(GR_Enhancers$peaks.signalValue, decreasing=TRUE)] 
    GR_Enhancers=GR_Enhancers[order(GR_Enhancers$peaks.qValue, decreasing=TRUE)] 
    
    regions_wider=resize(GR_Enhancers, width=2000, fix="center" )
    
    ##################Remove first ENCODE_blacklists from the data##############################################
    
    
   #remove blacklist regions
   test=as.matrix(findOverlaps(regions_wider, ENCODE_blacklist)) #query, subject
   if(nrow(test)!=0){
       GR_Enhancers=GR_Enhancers[-unique(test[,1])]
    }
    
    GR_Enhancers
    
}


distToTss <-function(peak, tss){
    
    peak <- resize(peak, width=1, fix="center")
    idx <- nearest(peak, tss) #(x, subject) returns indexes of the nearest subject from x
    naind=which(is.na(idx))
    notnaind=which(!is.na(idx))
    if(length(naind)!=0){
        idx=idx[-naind]
    }
    sgn <- as.numeric(ifelse(strand(tss)[idx] == "+", 1, -1))
    distances=data.frame(rep(5000000, length(peak)))
    distances[notnaind,]=(as.data.frame(ranges(peak))$start[notnaind] - as.data.frame(ranges(tss))$start[idx]) * sgn
}


ENCODE_blaclist_regions<-function(path_to_dir){
  
    # regions in the human genome that have anomalous, unstructured, high signal/read counts in NGS experiments independent of cell line and type of experiment
  
    #empirical from many ENCODE cell lines (not RNA-seq)
  
    #ultra-high signal artifact regions, 226
  #load blaclist regions
  path=paste(path_to_dir,"/Data/blacklists/",sep="")
  
  
  DAC<-read.table(paste(path,"wgEncodeDacMapabilityConsensusExcludable.bed.gz",sep=""), stringsAsFactors=FALSE) #411
  names(DAC)=c("chr", "start", "end", "name", "score", "strand")
  #[1] "High_Mappability_island" "Satellite_repeat"       
  #[3] "BSR/Beta"                "Low_mappability_island" 
  #[5] "(CATTC)n"                "LSU-rRNA_Hsa"           
  #[7] "centromeric_repeat"      "ALR/Alpha"              
  #[9] "SSU-rRNA_Hsa"            "telomeric_repeat"       
  #[11] "snRNA"                   "TAR1"                   
  #[13] "ACRO1"                   "chrM"                   
  #[15] "(GAGTG)n"                "(GAATG)n"     
  
  Duke<-read.table(paste(path,"wgEncodeDukeMapabilityRegionsExcludable.bed.gz",sep=""), stringsAsFactors=FALSE) #1649
  names(Duke)=c("chr", "start", "end", "name", "score", "strand")
  #[1] "TAR1"         "chrM"         "(GAATG)n"     "BSR/Beta"     "LSU-rRNA_Hsa"
  #[6] "(CATTC)n"     "(GAGTG)n"     "SSU-rRNA_Hsa" "ALR/Alpha"    "ACRO1"       
  #[11] "HSATII"      
  
  strandinformation=Rle( strand( rep( '+',nrow(DAC) ) ))
  DAC_GRanges=GRanges(seqnames = Rle(DAC$chr, rep(1, length(DAC$chr)) ), ranges = IRanges(start=DAC$start+1, end=DAC$end),  strand = strandinformation )
  
  strandinformation=Rle( strand( rep( '+',nrow(Duke) ) ))
  Duke_GRanges=GRanges(seqnames = Rle(Duke$chr, rep(1, length(Duke$chr)) ), ranges = IRanges(start=Duke$start+1, end=Duke$end),  strand = strandinformation )
  
  pathological=union(DAC_GRanges, Duke_GRanges) #1378
  pathological
  
}
  





explode_roi<-function(roi2, window, bin){
    
    GRanges(seqnames = Rle(as.character(seqnames(roi2[1])), rep(1, length(window/bin)) ), 
            ranges=IRanges(start=seq(start(roi2[1]), end(roi2[1])-bin+1, by=bin), 
                           end=seq(start(roi2[1])+bin-1, end(roi2[1]), bin) ), strand=Rle( strand( rep( '+',window/bin ) )))
}


get_coverage_from_bam<-function(roi, gal){
    assays(summarizeOverlaps(roi, gal, inter.feature=FALSE))$counts
    
    
}

round_correct<-function(value){
    
    #value can be a single value, vector or matrix
    #positive and negative values need to be handled separately
    
        
    pos_ind=which(value>=0)
    neg_ind=which(value<0)
        
    value[pos_ind]=trunc(value[pos_ind]+0.5)
    value[neg_ind]=trunc(value[neg_ind]-0.5)
    value
        
}

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





intervals_along_whole_chr<-function(chr_length, bin_size){
    bin_number=floor( chr_length/bin_size)
    starts<-1+(seq(1,bin_number,1)-1)*bin_size
    ends=1+( seq(1,bin_number,1) )*bin_size-1
    ends[bin_number]=chr_length
  
    
    
    strandinformation<-rep("+",bin_number)
    GRanges(seqnames = Rle(chr,  bin_number ), ranges = IRanges(start=starts, end = ends), strand = Rle(strand(strandinformation), rep(1,bin_number) ) )
    
}


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

enhancer_grouping_choose_type<-function(predictions_full, filename, type, threshold, enhancer_separation, overlap){
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




sample_random_locations<-function(chr,human.chromlens,widths){
  
  random_location=data.frame(chrom=chr, start=rep(0, length(widths)), end=rep(0, length(widths)))
  
  for(chrname in c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX")){
    #print(chrname)
    chr_ind=which(chr==chrname)
    location<-sample(1:(human.chromlens[chrname]),length(chr_ind), replace=FALSE)
    ind=which(((location-round(widths[chr_ind]/2))<0) || ((location+round(widths[chr_ind]/2))>human.chromlens[chrname]) )
    while(length(ind)!=0){
      location[ind]=sample(1:(human.chromlens[chrname]),length(ind), replace=FALSE)
      ind=which(((location-round(widths[chr_ind]/2))<0) || ((location+round(widths[chr_ind]/2))>human.chromlens[chrname]) )
    }
    
    random_location$start[chr_ind]=location-round(widths[chr_ind]/2)
    random_location$end[chr_ind]=location+(round(widths[chr_ind]/2)-1)
  }
  
  
  random_location
  
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



