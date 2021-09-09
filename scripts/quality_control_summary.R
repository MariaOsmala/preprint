library(optparse)

option_list = list(

make_option(c("-d", "--datadir"), type="character", default="",
            help="path to main folder [default= %default]", metavar="character")
)

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

path_to_dir=opt$datadir


cell_lines=c("K562", "GM12878")


#format:Filename<tab>numReads<tab>estFragLen<tab>corr_estFragLen<tab>PhantomPeak<tab>corr_phantomPeak<tab>argmin_corr<tab>min_corr<tab>Normalized SCC (NSC)<tab>Relative SCC (RSC)<tab>QualityTag)

#rm(data)
data<-data.frame()

# Below we cycle through directories. Remember the original working directory.
wd = getwd()
for(cl in cell_lines){

  path=paste0(path_to_dir,"/",cl,"/phantompeakqualtools/")

  # The path may not exist if the cell line is not being processed.
  # Only try reading if the folder exists.
  if(file.exists(path)){
    setwd(path)

    mods=dir(pattern=".out")

    for(mod in mods){

      tmp=try(read.table(mod, stringsAsFactors=FALSE),silent=TRUE)
      if(class(tmp)!="try-error"){

        data=rbind(data, cbind(cl,tmp))

      }
    }

    # Go back to the original working directory
    setwd(wd)
  }
}

colnames(data)=c("cell_line","name", "#mapped_reads", "estFragLen", "corr_estFragLen", "phantomPeak", "corr_phantomPeak", "argmin_corr", "min_corr", "NSC=corr_estFraglen/min_corr >=1.05", "RSC=(corr_estFraglen - min_corr)/(corr_phantomPeak - min_corr) >=0.8", "QualityTag based on RSC -2:veryLow,-1:Low,0:Medium,1:High,2:veryHigh")

fragLengths=lapply(list(data$estFragLen), FUN=strsplit, split=",")
corrs=lapply(list(data$corr_estFragLen), FUN=strsplit, split=",")

how_many=max(unlist(lapply(fragLengths[[1]], length)))


for(i in 1:length(fragLengths[[1]])){

  if(length(fragLengths[[1]][[i]])< how_many){

    fragLengths[[1]][[i]]=as.numeric( c(fragLengths[[1]][[i]], rep("NA", how_many - length(fragLengths[[1]][[i]]))  ) )
    corrs[[1]][[i]]=as.numeric( c(corrs[[1]][[i]], rep("NA", how_many - length(corrs[[1]][[i]]))  ) )

  }

  fragLengths[[1]][[i]]=as.numeric(fragLengths[[1]][[i]])
  corrs[[1]][[i]]=as.numeric(corrs[[1]][[i]])

}
new_data<-cbind(data[,1:3], do.call(rbind, fragLengths[[1]]), do.call(rbind, corrs[[1]]), data[, 6:12])
colnames(new_data)=c("cell_line","name", "#mapped_reads", "1st_estFragLen","2nd_estFragLen","3rd_estFragLen", "1st_corr_estFragLen","2nd_corr_estFragLen","3rd_corr_estFragLen", "phantomPeak", "corr_phantomPeak", "argmin_corr", "min_corr", "NSC=corr_estFraglen/min_corr >=1.05", "RSC=(corr_estFraglen - min_corr)/(corr_phantomPeak - min_corr) >=0.8", "QualityTag based on RSC -2:veryLow,-1:Low,0:Medium,1:High,2:veryHigh")



max_value=apply(new_data[,7:9],1,max,na.rm=TRUE)

max_index=apply(new_data[,7:9],1,which.max) #all are ones


write.table(new_data, file=paste(path_to_dir,"/phantompeakqualtools.txt",sep=""), quote=FALSE, row.names=FALSE, sep="\t")
