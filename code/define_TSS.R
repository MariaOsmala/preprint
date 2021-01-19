#library(Rsamtools)
#library(snow)
#library(spp)
#library(accelerometry)
#library(Biostrings)
#library(bitops)
#library(BSgenome.Hsapiens.UCSC.hg19)
#library(circlize)
#library(doParallel)
#library(foreach)
#library(gdata)
#library(GenomicRanges)
#library(GetoptLong)
#library(ggplot2)
#library(grid)
#library(gridExtra)
#library(MASS)
library(optparse)
#library(pryr)
#library(RColorBrewer)
#library(reshape2)
#library(ROCR)
#library(rtracklayer)
#library(ShortRead)
#library(stringr)


option_list = list(
  
  make_option("--pathToDir", type="character", default="", 
              help="path to main folder [default= %default]", metavar="character")
); 


opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

path=opt$pathToDir

Gencode_TSS=read.table(paste(path,"/GENCODE_TSS/gencode.v27lift37.annotation.gtf.gz",sep=""), header=FALSE, 
                       sep="\t", stringsAsFactors=FALSE) # 2629275 x 9

names(Gencode_TSS)=c("chrom","annotation","feature_type","start","end","score","strand","genomic_phase","info")

extract_info_attributes <- function(info){
  
  gsub(";","",do.call(rbind, strsplit(unlist(strsplit(info[1],"; "))," ")))[,1]
}

attribute_names=unique(unlist(lapply(as.list(Gencode_TSS$info), extract_info_attributes))) #24


list_to_matrix <- function(list){
  
  do.call(rbind, list)
}

attribute_data=lapply(lapply(lapply(lapply(lapply(as.list(Gencode_TSS$info),strsplit, "; "), unlist), strsplit, " "), 
                             list_to_matrix),gsub,pattern=";",replacement="")

find_attribute_value <-function(table, atr){
  
  table[which(table[,1]==atr),2]
  
}

#remove tag from attribute_names, handle it separately if needed, tags contain CAGE_supported_TSS


#there can be 8 different tags

atr="tag"
tmp=lapply(attribute_data, find_attribute_value, atr=atr)


extend_vector <- function(x, n, fill_value=""){
  xx <- rep(x, length.out=n)
  if(length(x) < n){
    xx[(length(x)+1):n] <- fill_value 
  }
  xx
}

tmp3=do.call(rbind, lapply(lapply(lapply(tmp, strsplit, " "),unlist), extend_vector, n=8))
tmp4=as.data.frame(tmp3, stringsAsFactors=FALSE)
names(tmp4)=c("tag_1","tag_2","tag_3","tag_4","tag_5","tag_6","tag_7","tag_8")


Gencode_TSS=cbind(Gencode_TSS, tmp4)

attribute_names=attribute_names[-which(attribute_names=="tag")]
attribute_names=attribute_names[-which(attribute_names=="ont")]

for(atr in attribute_names){
  print(atr)  
  tmp=lapply(attribute_data, find_attribute_value, atr=atr)
  nonempty_ind=which(unlist(lapply(tmp, length))==1)
  Gencode_TSS[, atr]=""
  Gencode_TSS[nonempty_ind, atr]=unlist(tmp)
  
}

Gencode_TSS=Gencode_TSS[,-9] #remove the info

save.image(paste(path,"GENCODE_TSS/GENCODE.RData",sep=""))


#chrom
allowed_chroms=c("chr1","chr2","chr3","chr4","chr5","chr6", "chr7","chr8","chr9","chr10",     
                 "chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX","chrY","chrM")
dd=Gencode_TSS$chrom %in% allowed_chroms
Gencode_TSS=Gencode_TSS[which(dd==TRUE),]


#[3] "feature_type": "gene"           "transcript"     "exon"           "CDS"           
# "start_codon"    "stop_codon"     "UTR"            "Selenocysteine"

#"exon_number" #all except gene/transcript/Selenocysteine	integer (exon position in the transcript from its 5' end)                     
#"exon_id"                          

#consider only transcripts

Gencode_transcripts=Gencode_TSS[which(Gencode_TSS$feature_type=="transcript"),]

#GRanges object with TSS strand information present
#add 1 to the start coordinates so we can get the 1-based coordinate system

#GFF are 1 based coordinate system

#returns TSS coordinates in strand-specific and strand-nonspecific format
annotation_to_GRanges <- function(annotation){
  GR=GRanges( seqnames = Rle(annotation$chrom, rep(1, nrow(annotation) ) ), 
              ranges = IRanges(start=annotation$start, end=annotation$end), 
              strand = Rle( strand( annotation$strand ), rep(1,nrow(annotation) ) )) 
  
  #TSS as a single nucleotide, resize takes the strand automatically into account
  GR=resize(GR,1,fix="start")
  GR=unique(GR)
  
  #TSS defined in positive strand only
  GR_positive=GR
  strand(GR_positive)=strand( rep("+", length(GR_positive) ))
  
  result <- list(GR=GR,GR_positive=GR_positive)
  return(result)
  
}


result=annotation_to_GRanges(Gencode_transcripts)

GR_Gencode_TSS=result$GR
GR_Gencode_TSS_positive=result$GR_positive

#define also Gencode_protein_coding_TSS
Gencode_protein_coding=Gencode_transcripts[which(Gencode_transcripts$transcript_type=="protein_coding"),] #81326
result=annotation_to_GRanges(Gencode_protein_coding)
GR_Gencode_protein_coding_TSS=result$GR
GR_Gencode_protein_coding_TSS_positive=result$GR_positive
rm(result)

saveRDS(GR_Gencode_protein_coding_TSS, paste(path,"/GENCODE_TSS/","GR_Gencode_protein_coding_TSS.RDS",sep=""))
saveRDS(GR_Gencode_protein_coding_TSS_positive,paste(path,"/GENCODE_TSS/","GR_Gencode_protein_coding_TSS_positive.RDS",sep="")) #This is used in the subsequent analysis
saveRDS(GR_Gencode_TSS,paste(path,"/GENCODE_TSS/","GR_Gencode_TSS.RDS",sep=""))
saveRDS(GR_Gencode_TSS_positive,paste(path,"/GENCODE_TSS/","GR_Gencode_TSS_positive.RDS",sep=""))
