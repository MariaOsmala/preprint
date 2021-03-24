library(GenomicRanges)

find_random_with_signal <- function(coverage, regions, threshold = 5, window = 1000, N = NULL,
                                    p300 = NULL, max_dist_to_p300 = 2000,
                                    TSS = NULL, max_dist_to_TSS = 2000,
                                    blacklist = NULL)
{
    print('Loading data...')
    # Select the $data parts of the coverage list
    coverage <- unlist(coverage, recursive = FALSE)[c(TRUE, FALSE)]

    # Concatenate to one big matrix
    coverage <- t(do.call(rbind, coverage))

    # Don't count negative coverage values 
    coverage[coverage < 0] <- 0

    # Select all regions which have enough coverage
    accepted <- rowSums(coverage) >= threshold
    regions <- regions[accepted]
    print(paste0('Selected ', length(regions), ' regions with enough coverage.'))

    # Blacklist to remove problematic regions
    if (is.null(blacklist)) {
        blacklist <- GRanges()
    }

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

    # Remove all blacklisted regions
    to_remove <- unique(from(findOverlaps(regions, blacklist)))
    if (length(to_remove) > 0) {
        regions <- regions[-to_remove]
    }
    print(paste0("#regions after blacklist: ", length(regions)))

    # Select N random locations, max 1 location per region
    if (is.null(N)) {
        N <- length(regions)
    }
    regions <- sort(sample(regions, N))
    samples <- sample(seq(1, bin_size), N, replace = TRUE)
    start(regions) <- start(regions) + samples
    width(regions) <- 1

    regions
}
