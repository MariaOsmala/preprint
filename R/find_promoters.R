#' Find promoter sites suitable to serve as training data.
#' @importFrom GenomicRanges mcols
#' @importFrom IRanges from to resize findOverlaps

#' Find promoter sites
#'
#' @description
#' Finds promoter sites suitable to serve as training data.
#'
#' @details
#' Promoter sites are sites on the genome that have TSS annotations and are
#' also DNase binding sites.
#'
#' @param TSS_annotations
#'     A `GRanges` object containing the TSS annotations that will be used to
#'     mark candidate promoter sites.
#' @param DNase 
#'     A `GRanges` object containing DNase binding sites. There should be a
#'     `peak` metadata column here, At least three metadata columns should be
#'     present: `peak`, `qValue` and `signalValue`.
#' @param window
#'     The length of the window to cut around each promoter site. This window
#'     can later be used to create the corresponding coverage profile using
#'     `create_profiles`. The window will be centered on the promoter site.
#' @param N
#'     The number of promoter sites to find. Promoter sites with the highest
#'     qValue will be used (using signalValue to resolve ties).
#' @param between_promoter_distance
#'     Desired minimum distance between the selected promoter sites.
#' @param blacklist
#'     An optional `GRanges` object containing areas that should be avoided
#'     when searching for suitable promoter sites.
#' @param verbose
#'     A logical scalar indicating whether to print out informational messages
#'     during the search for suitable promoter sites.
#'
#' @return
#'     A `GRanges` object containing windows centered around suitable promoter
#'     sites. The promoter site itself lies in the center of the window. The
#'     `RGanges` object has a `type` metadata column that is set to
#'     `'promoter'`.
#'
#' @seealso [find_enhancers()] for finding suitable enhancer sites and
#'          [find_random()] for finding suitable random sites.
#' @export
find_promoters <- function(TSS_annotations, DNase, window = 1000, N = NULL,
                           between_promoter_distance = 2000, blacklist = NULL,
                           verbose = TRUE)
{
    if (verbose) cat('Finding suitable promoter sites:\n')

    # Find the closest overlap between DNase and TSS_annotations
    dist <- GenomicRanges::distanceToNearest(DNase, TSS_annotations, ignore.strand = TRUE)
    TSS_with_DNase <- TSS_annotations[to(dist)]
    mcols(TSS_with_DNase) <- mcols(DNase[from(dist)])[c('signalValue', 'pValue')]
    mcols(TSS_with_DNase)$distance <- mcols(dist)$distance

    # Compute the distance from the peak in the DNase to the corresponding TSS_annotations
    DNase_peak <- DNase
    start(DNase_peak) <- start(DNase_peak) + DNase_peak$peak
    width(DNase_peak) <- 1
    dist_peak <- GenomicRanges::distanceToNearest(DNase_peak, TSS_annotations, ignore.strand = TRUE)
    mcols(TSS_with_DNase)$peakDistance <- mcols(dist_peak[from(dist)])$distance

    # Promoter sites are TSS that overlap with DNase
    TSS_with_DNase <- TSS_with_DNase[order(mcols(TSS_with_DNase)$distance, decreasing = FALSE)]
    to_keep <- from(findOverlaps(TSS_with_DNase, DNase, ignore.strand = TRUE))
    promoters <- TSS_with_DNase[to_keep] #60528

    promoters <- promoters[order(promoters$distance)]
    tmp <- resize(promoters, between_promoter_distance / 2, fix = 'center')
    tmp <- findOverlaps(tmp, ignore.strand = TRUE)
    tmp <- tmp[from(tmp) != to(tmp)]  # Ranges always overlap themselves
    to_drop <- unique(c(from(tmp), to(tmp)))  # All ranges that overlap with any other ranges
    if (length(to_drop) > 0) {
        promoters <- promoters[-to_drop]
    }
    
    if (verbose) cat(paste0('    #promoters that are nicely isolated: ', length(promoters), '\n'))

    # TODO: not sure why promoter sites that didn't overlap with TSS_annotations were kept until now
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
