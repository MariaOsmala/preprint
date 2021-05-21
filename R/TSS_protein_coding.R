#' @importFrom IRanges resize
#' @importFrom GenomicRanges mcols strand<-

#' Use the histone reads to create profiles for a list of genomic sites
#' @export
TSS_protein_coding <- function(gencode_tss_file, verbose = TRUE)
{
    if (verbose) cat('Reading TSS annotations...')
    gencode <- rtracklayer::import(gencode_tss_file)
    gencode <- resize(gencode, 1, fix = "start")  # We only care about the start of the read

    # Make selections
    chroms_of_interest <- c(sapply(1:22, function(x) paste0('chr', x)), "chrX", "chrY", "chrM")
    gencode <- GenomeInfoDb::keepSeqlevels(gencode, chroms_of_interest, pruning.mode = "coarse")
    gencode_transcripts <- gencode[mcols(gencode)$type == "transcript"]
    gencode_protein_coding <- gencode_transcripts[mcols(gencode_transcripts)$transcript_type == "protein_coding"]

    # As we no longer have range lengths, there are duplicates
    gencode_protein_coding <- unique(gencode_protein_coding)

    # Strip out metadata
    mcols(gencode_protein_coding) <- NULL

    # Ignore strands
    strand(gencode_protein_coding) <- "*"

    if (verbose) cat(' done.\n')

    gencode_protein_coding
}
