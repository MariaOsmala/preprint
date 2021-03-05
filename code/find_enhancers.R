library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg19)

find_enhancers  <- function(p300, DNase, window = 1000, N = NULL,
                            promoters = NULL, max_dist_to_promoter = 2000,
                            blacklist = NULL)
{
    strand(p300) <- "*"  # This is important

    # Trim the p300 ranges to just the peak values
    start(p300) <- start(p300) + p300$peak
    width(p300) <- 1

    print(paste0("#p300 peaks: ", length(p300)))

     # TODO: why does chrM needs to be removed?
    p300 <- p300[seqnames(p300) != "chrM"]

    if (!is.null(promoters)) { 
        # Remove p300 peaks that are too close to a promotor range, as specified by
        # max_dist_to_promoter.
        # We ignore the strand, which means all strands are presumed to be "+".
        dist <- distanceToNearest(p300, promoters, ignore.strand = TRUE)
        to_remove <- from(dist)[mcols(dist)$distance < (max_dist_to_promoter - 1)]  # FIXME: why the -1?
        p300 <- p300[-to_remove]
    }
    print(paste0("#p300 peaks away from promoters: ", length(p300)))

    # Enhancers are p300 peaks that overlap with DNase ranges
    enhancers <- p300[unique(from(findOverlaps(p300, DNase)))]
    print(paste0("#enhancers: ", length(enhancers)))

    # Cut a window around the peak.
    # If window%%2==0, add 1 to window, extend the enhancer window/2 downstream and window/2 upsteam
    # If window%%2==1, extend the enhancer window (window-1)/2 downstream and upstream
    enhancers <- resize(enhancers, width = window + (window+1)%%2, fix="center")

    # Remove enhancers which window falls outside the bounds of the sequence
    human.chromlens <- seqlengths(Hsapiens)
    seqlengths(enhancers) <- human.chromlens[seqnames(seqinfo(enhancers))]
    seq_lengths <- as.numeric(seqlengths(enhancers)[as.vector(seqnames(enhancers))])
    enhancers <- enhancers[start(enhancers) >= 1 & end(enhancers) <= seq_lengths]

    if (!is.null(blacklist)) {
        # Remove enhancers that fall within blacklisted regions
        enhancers <- enhancers[-unique(from(findOverlaps(enhancers, blacklist)))]
    }
    print(paste0("#enhancers after blacklist: ", length(enhancers)))

    if (!is.null(N)) {
        # Sort first by qValue, then by signalValue
        selection <- order(enhancers$qValue, enhancers$signalValue, decreasing = TRUE)[1:N]
        enhancers <- enhancers[selection]
    }
    print(paste0("#enhancers after selecting N: ", length(enhancers)))

    enhancers
}
