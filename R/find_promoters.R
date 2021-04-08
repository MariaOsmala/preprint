#' Find promoter sites suitable to serve as training data.
#' @importFrom GenomicRanges mcols
#' @importFrom IRanges from to resize findOverlaps
#' @export
find_promoters <- function(TSS_annotation, DNase, window = 1000, N = NULL,
                           between_TSS_distance = 2000, blacklist = NULL,
                           verbose = TRUE)
{
    if (verbose) cat('Finding suitable promoter sites:\n')

    # Find the closest overlap between DNase and TSS_annotation
    dist <- GenomicRanges::distanceToNearest(DNase, TSS_annotation, ignore.strand = TRUE)
    TSS_with_DNase <- TSS_annotation[to(dist)]
    mcols(TSS_with_DNase) <- mcols(DNase[from(dist)])[c('signalValue', 'pValue')]
    mcols(TSS_with_DNase)$distance <- mcols(dist)$distance

    # Compute the distance from the peak in the DNase to the corresponding TSS_annotation
    DNase_peak <- DNase
    start(DNase_peak) <- start(DNase_peak) + DNase_peak$peak
    width(DNase_peak) <- 1
    dist_peak <- GenomicRanges::distanceToNearest(DNase_peak, TSS_annotation, ignore.strand = TRUE)
    mcols(TSS_with_DNase)$peakDistance <- mcols(dist_peak[from(dist)])$distance

    # Order
    TSS_with_DNase <- TSS_with_DNase[order(mcols(TSS_with_DNase)$distance, decreasing = FALSE)]

    # Promoter sites are TSS that overlap with DNase
    to_keep <- from(findOverlaps(TSS_with_DNase, DNase, ignore.strand = TRUE))
    promoters <- TSS_with_DNase[to_keep] #60528

    promoters <- promoters[order(promoters$distance)]
    tmp <- resize(promoters, between_TSS_distance / 2, fix = 'center')
    tmp <- findOverlaps(tmp, ignore.strand = TRUE)
    tmp <- tmp[from(tmp) != to(tmp)]  # Ranges always overlap themselves
    to_drop <- unique(c(from(tmp), to(tmp)))  # All ranges that overlap with any other ranges
    if (length(to_drop) > 0) {
        promoters <- promoters[-to_drop]
    }
    
    if (verbose) cat(paste0('    #promoters that are nicely isolated: ', length(promoters), '\n'))

    # ???
    promoters <- promoters[promoters$distance == 0]

    # Create a window around the promoters
    promoters <- resize(promoters, width = window, fix = 'center' )

    if (!is.null(blacklist)) {
        # Remove promoters that fall within blacklisted regions
        to_remove <- unique(from(findOverlaps(promoters, blacklist)))
        if (length(to_remove) > 0) {
            promoters <- promoters[-to_remove]
        }
        
        if (verbose) cat(paste0('    #promoters after blacklist: ', length(promoters), '\n'))
    }

    if (!is.null(N)) {
        # Sort first by qValue, then by signalValue, then by peakDistance
        promoters <- promoters[order(promoters$pValue, promoters$signalValue, decreasing = TRUE)]
        promoters <- promoters[order(promoters$peakDistance, decreasing = FALSE)]
        promoters <- promoters[1:N]
        if (verbose) cat(paste0('    #promoters after selecting N: ', length(promoters), '\n'))
    }

    mcols(promoters)$type <- as.factor('promoter')
    promoters
}
