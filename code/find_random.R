library(BSgenome.Hsapiens.UCSC.hg19)
library(GenomicRanges)

find_random <- function(window = 1000, N = NULL,
                        p300 = NULL, max_dist_to_p300 = 2000,
                        TSS = NULL, max_dist_to_TSS = 2000,
                        blacklist = NULL)
{
    allowed_chroms=c("chr1",  "chr2",  "chr3",  "chr4",  "chr5",  "chr6",  "chr7",  "chr8",  "chr9", 
    "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18",
    "chr19", "chr20", "chr21", "chr22", "chrX") 
    whole_genome <- GRanges(allowed_chroms, IRanges(start=1, end=seqlengths(Hsapiens)[allowed_chroms]))

    if (is.null(blacklist)) {
        blacklist = GRanges()
    }

    # Remove regions that are too close to a p300 peak, as specified by
    # max_dist_to_p300
    if (!is.null(p300)) {
        p300_exclusion_zone <- p300
        p300_exclusion_zone <- resize(p300_exclusion_zone, width = max_dist_to_p300 + (max_dist_to_p300 + 1) %% 2, fix="center")
        blacklist <- union(blacklist, p300_exclusion_zone)
    }

    # Remove regions that are too close to an TSS annotation, as specified by
    # max_dist_to_TSS
    if (!is.null(TSS)) {
        TSS_exclusion_zone <- TSS
        TSS_exclusion_zone <- resize(TSS_exclusion_zone, width = max_dist_to_TSS + (max_dist_to_TSS + 1) %% 2, fix="center")
        blacklist <- union(blacklist, TSS_exclusion_zone)
    }

    blacklist <- keepSeqlevels(blacklist, allowed_chroms, pruning.mode = "coarse")

    # reduce() will greatly reduce the number of ranges by merging adjacent ranges.
    blacklist <- reduce(blacklist)

    print("Finding random regions that satisfy the criteria.")
    random_regions <- regioneR::createRandomRegions(N, length.mean = window, length.sd = 0, genome = whole_genome, mask = blacklist, non.overlapping = FALSE)
    random_regions <- sort(random_regions)
    mcols(random_regions)$type <- as.factor('random')

    random_regions
}
