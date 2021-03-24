library(GenomicRanges)

find_random_with_signal <- function(coverage, regions, threshold = 5, window = 1000, N = NULL,
                                    p300 = NULL, max_dist_to_p300 = 2000,
                                    TSS = NULL, max_dist_to_TSS = 2000,
                                    blacklist = NULL)
{
    # All blacklisted regions will be omitted
    if (is.null(blacklist)) {
        blacklist <- GRanges()
    }

    # Select the $data parts of the coverage list
    coverage <- unlist(coverage, recursive = FALSE)[c(TRUE, FALSE)]

    # Concatenate to one big matrix
    coverage <- t(do.call(rbind, coverage))

    # Don't count negative coverage values 
    coverage[coverage < 0] <- 0

    # Blacklist all regions which don't have enough coverage
    blacklist <- union(blacklist, regions[rowSums(coverage) < threshold])

    # Remove regions that are too close to a p300 peak, as specified by
    # max_dist_to_p300
    if (!is.null(p300)) {
        dist <- distanceToNearest(regions, p300, ignore.strand = TRUE)
        too_close <- from(dist)[mcols(dist)$distance < (max_dist_to_p300 - 1)]  # FIXME: why the -1?
        blacklist <- union(blacklist, regions[too_close])
    }

    # Remove regions that are too close to an TSS annotation, as specified by
    # max_dist_to_TSS
    if (!is.null(TSS)) {
        dist <- distanceToNearest(regions, TSS, ignore.strand = TRUE)
        too_close <- from(dist)[mcols(dist)$distance < (max_dist_to_TSS - 1)]  # FIXME: why the -1?
        blacklist <- union(blacklist, regions[too_close])
    }

    allowed_chroms=c("chr1",  "chr2",  "chr3",  "chr4",  "chr5",  "chr6",  "chr7",  "chr8",  "chr9", 
    "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18",
    "chr19", "chr20", "chr21", "chr22", "chrX") 
    blacklist <- keepSeqlevels(blacklist, allowed_chroms, pruning.mode = "coarse")

    regioneR::createRandomRegions(N, length.mean = window, length.sd = 0, genome = regions, mask = blacklist)
}
