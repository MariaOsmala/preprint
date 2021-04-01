library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg19)

find_enhancers  <- function(p300, DNase, window = 1000, N = NULL,
                            TSS = NULL, max_dist_to_promoter = 2000,
                            blacklist = NULL)
{
    strand(p300) <- "*"  # This is important

    # Trim the p300 ranges to just the peak values
    start(p300) <- start(p300) + p300$peak
    width(p300) <- 1

    print(paste0("#p300 peaks: ", length(p300)))

     # TODO: why does chrM needs to be removed?
    p300 <- p300[seqnames(p300) != "chrM"]

    if (!is.null(TSS)) { 
        # Remove p300 peaks that are too close to a promotor range, as specified by
        # max_dist_to_promoter.
        # We ignore the strand, which means all strands are presumed to be "+".
        dist <- distanceToNearest(p300, TSS, ignore.strand = TRUE)
        to_drop <- from(dist)[mcols(dist)$distance < (max_dist_to_promoter - 1)]  # FIXME: why the -1?
        if (length(to_drop) > 0) {
            p300 <- p300[-to_drop]
        }
    }
    print(paste0("#p300 peaks away from promoters: ", length(p300)))

    # Enhancers are p300 peaks that overlap with DNase ranges
    enhancers <- p300[unique(from(findOverlaps(p300, DNase)))]
    print(paste0("#enhancers: ", length(enhancers)))

    # Cut a window around the peak
    enhancers <- resize(enhancers, width = window, fix="center")

    # Remove enhancers which window falls outside the bounds of the sequence
    human.chromlens <- seqlengths(Hsapiens)
    seqlengths(enhancers) <- human.chromlens[seqnames(seqinfo(enhancers))]
    seq_lengths <- as.numeric(seqlengths(enhancers)[as.vector(seqnames(enhancers))])
    enhancers <- enhancers[start(enhancers) >= 1 & end(enhancers) <= seq_lengths]

    if (!is.null(blacklist)) {
        # Remove enhancers that fall within blacklisted regions
        to_drop <- unique(from(findOverlaps(enhancers, blacklist)))
        if (length(to_drop) > 0) {
            enhancers <- enhancers[-to_drop]
        }
        print(paste0("#enhancers after blacklist: ", length(enhancers)))
    }

    if (!is.null(N)) {
        # Sort first by qValue, then by signalValue
        selection <- order(enhancers$qValue, enhancers$signalValue, decreasing = TRUE)[1:N]
        enhancers <- enhancers[selection]
    }
    print(paste0("#enhancers after selecting N: ", length(enhancers)))

    mcols(enhancers)$type <- as.factor('enhancer')
    enhancers
}
