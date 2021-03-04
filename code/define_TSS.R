library(GenomicRanges)
library(optparse)

option_list = list(
  make_option("--pathToDir", type="character", default="",
              help="path to main folder [default= %default]", metavar="character")
);
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

path=opt$pathToDir

gencode <- rtracklayer::import(paste0(path, "/GENCODE_TSS/gencode.v27lift37.annotation.gtf.gz"))
gencode <- resize(gencode, 1, fix = "start")  # We only care about the start of the read

# Make selections
allowed_chroms <- c(sapply(1:22, function(x) paste0('chr', x)), "chrX", "chrY", "chrM")
gencode <- gencode[chrom(gencode) %in% allowed_chroms]
gencode_transcripts <- gencode[mcols(gencode)$type == "transcript"]
gencode_protein_coding <- gencode_transcripts[mcols(gencode_transcripts)$transcript_type == "protein_coding"]

# As we no longer have range lengths, there are duplicates
gencode_transcripts <- unique(gencode_transcripts)  
gencode_protein_coding <- unique(gencode_protein_coding)

# Strip out metadata
mcols(gencode_transcripts) <- NULL
mcols(gencode_protein_coding) <- NULL

# Make versions where all strands are positive
gencode_transcripts_positive <- gencode_transcripts
gencode_protein_coding_positive <- gencode_protein_coding
strand(gencode_transcripts_positive) <- "*"
strand(gencode_protein_coding_positive) <- "*"

# Save everything
saveRDS(gencode_transcripts, paste0(path, "/GENCODE_TSS/","GR_Gencode_TSS2.RDS"))
saveRDS(gencode_transcripts_positive, paste0(path, "/GENCODE_TSS/","GR_Gencode_TSS_positive2.RDS"))
saveRDS(gencode_protein_coding, paste0(path, "/GENCODE_TSS/","GR_Gencode_protein_coding_TSS2.RDS"))
saveRDS(gencode_protein_coding_positive, paste0(path, "/GENCODE_TSS/","GR_Gencode_protein_coding_TSS_positive2.RDS"))
