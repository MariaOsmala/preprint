library(GenomicRanges)
library(optparse)

option_list = list(
  make_option("--pathToDir", type="character", default="",
              help="path to main folder [default= %default]", metavar="character")
);
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

path=opt$pathToDir

Gencode_TSS <- read.table(paste0(path, "/GENCODE_TSS/gencode.v27lift37.annotation.gtf.gz"),
                          header=FALSE,
                          sep="\t",
                          stringsAsFactors=FALSE) # 2629275 x 9

names(Gencode_TSS) <- c("chrom", "annotation", "feature_type", "start", "end",
                        "score", "strand", "genomic_phase", "info")

# Parsing the info column into a list of lists.
#
# The info column of Gencode_TSS is one big string with key-value pairs like this:
# "gene_id ENSG00000223972.5_2; gene_type transcribed_unprocessed_pseudogene; gene_name DDX11L1"
# We first split this string into words. The even words are the values, and the
# uneven words are the names for those values. A complication is that there are
# multiple values for "tag", which we rename "tag_1", "tag_2", etc.
parse_info <- function(info_str) {
  info_words <- strsplit(info_str, " ")[[1]]
  info_list <- info_words[seq(2, length(info_words), by=2)]
  info_list <- sub(";", "", info_list)
  info_names <- info_words[seq(1, length(info_words), by=2)]
  tags = which(info_names == "tag")
  for (i in seq_along(tags)) {
    info_names[tags[i]] <- paste0("tag", "_", i)
  }
  names(info_list) <- info_names
  return(info_list)
}

gen_info <- lapply(Gencode_TSS$info, parse_info)

# Turn gen_info from a list of lists into a data.frame. A complication here is
# that not all the lists define the same values. We must stratify them first by
# indexing them by all possible names, which will create N/A's for any missing
# values. We first pass the data through the matrix() function to create a
# data.frame where all values are of the same type (character).
attribute_names <- unique(unlist(lapply(gen_info, names)))
gen_info <- data.frame(
  matrix(
    unlist(lapply(gen_info, function(x) x[attribute_names])),
    nrow=length(gen_info),
    byrow=TRUE)
  )
names(gen_info) <- attribute_names

# We don't need this column
gen_info$ont <- NULL

Gencode_TSS <- cbind(Gencode_TSS, gen_info)
save.image(paste0(path, "/GENCODE_TSS/GENCODE.RData"))

allowed_chroms <- c(sapply(1:22, function(x) paste0('chr', x)),
                    "chrX", "chrY", "chrM")
dd <- Gencode_TSS$chrom %in% allowed_chroms
Gencode_TSS <- Gencode_TSS[dd,]


#[3] "feature_type": "gene"           "transcript"     "exon"           "CDS"
# "start_codon"    "stop_codon"     "UTR"            "Selenocysteine"

#"exon_number" #all except gene/transcript/Selenocysteine	integer (exon position in the transcript from its 5' end)
#"exon_id"

#consider only transcripts

Gencode_transcripts <- Gencode_TSS[Gencode_TSS$feature_type == "transcript",]

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
