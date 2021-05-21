#' @importFrom IRanges end findOverlaps from resize start start<- width width<-
#' @importFrom GenomicRanges strand<- seqnames seqinfo mcols mcols<-
#' @importFrom GenomeInfoDb seqlengths seqlengths<-

#' Find enhancer sites
#'
#' @description
#' Finds enhancer sites suitable to serve as training data.
#'
#' @details
#' Enhancer sites are sites on the genome that correspond to p300 peaks and are
#' also DNase binding sites.
#'
#' @param p300 
#'     A `GRanges` object containing the P300 data. At least three metadata
#'     columns should be present: `peak`, `qValue` and `signalValue`.
#' @param DNase 
#'     A `GRanges` object containing DNase binding sites.
#' @param window
#'     The length of the window to cut around each enhancer site. This window
#'     can later be used to create the corresponding coverage profile using
#'     `create_profiles`. The window will be centered on the enhancer site.
#' @param N
#'     The number of enhancer sites to find. Enhancer sites with the highest
#'     qValue will be used (using signalValue to resolve ties).
#' @param TSS_annotations
#'     When specified, enhancer sites that are too close to a promoter site
#'     will not be used. Set this to a `GRanges` object containing the TSS
#'     annotations.
#' @param min_dist_to_promoter
#'     When the `TSS` parameter is specified, this parameter controls the
#'     minimum distance that enhancer sites should have to promoter sites.
#' @param blacklist
#'     An optional `GRanges` object containing areas that should be avoided
#'     when searching for suitable enhancer sites.
#' @param verbose
#'     A logical scalar indicating whether to print out informational messages
#'     during the search for suitable enhancer sites.
#'
#' @return
#'     A `GRanges` object containing windows centered around suitable enhancer
#'     sites. The enhancer site itself lies in the center of the window. The
#'     `RGanges` object has a `type` metadata column that is set to
#'     `'enhancer'`.
#'
#' @seealso [find_promoters()] for finding suitable promoter sites and
#'          [find_random()] for finding suitable random sites.
#' @export
find_enhancers  <- function(p300, DNase, window = 1000, N = NULL,
                            TSS_annotations = NULL, min_dist_to_promoter = 2000,
                            blacklist = NULL, verbose = TRUE)
{
    if (verbose) cat('Finding suitable enhancer sites:\n')
    strand(p300) <- '*'  # This is important

    # Trim the p300 ranges to just the peak values
    start(p300) <- start(p300) + p300$peak
    width(p300) <- 1

    if (verbose) cat(paste0('    #p300 peaks: ', length(p300), '\n'))

    # TODO: why does chrM need to be removed?
    p300 <- p300[seqnames(p300) != 'chrM']

    if (!is.null(TSS_annotations)) { 
        # Remove p300 peaks that are too close to a promotor range, as specified by
        # min_dist_to_promoter.
        # We ignore the strand, which means all strands are presumed to be '+'.
        dist <- GenomicRanges::distanceToNearest(p300, TSS_annotations, ignore.strand = TRUE)
        to_drop <- from(dist)[mcols(dist)$distance < (min_dist_to_promoter - 1)]  # FIXME: why the -1?
        if (length(to_drop) > 0) {
            p300 <- p300[-to_drop]
        }
    }
    
    if (verbose) cat(paste0('    #p300 peaks away from promoters: ', length(p300), '\n'))

    # Enhancers are p300 peaks that overlap with DNase ranges
    enhancers <- p300[unique(from(findOverlaps(p300, DNase)))]
    if (verbose) cat(paste0('    #enhancers found: ', length(enhancers), '\n'))

    # Cut a window around the peak
    enhancers <- resize(enhancers, width = window, fix='center')

    # Remove enhancers which window falls outside the bounds of the sequence
    human.chromlens <- seqlengths(BSgenome.Hsapiens.UCSC.hg19::Hsapiens)
    seqlengths(enhancers) <- human.chromlens[seqnames(seqinfo(enhancers))]
    seq_lengths <- as.numeric(seqlengths(enhancers)[as.vector(seqnames(enhancers))])
    enhancers <- enhancers[start(enhancers) >= 1 & end(enhancers) <= seq_lengths]

    if (!is.null(blacklist)) {
        # Remove enhancers that fall within blacklisted regions
        to_drop <- unique(from(findOverlaps(enhancers, blacklist)))
        if (length(to_drop) > 0) {
            enhancers <- enhancers[-to_drop]
        }
        if (verbose) cat(paste0('    #enhancers after blacklist: ', length(enhancers), '\n'))
    }

    if (!is.null(N)) {
        # Sort first by qValue, then by signalValue
        selection <- order(enhancers$qValue, enhancers$signalValue, decreasing = TRUE)[1:N]
        enhancers <- enhancers[selection]
        if (verbose) cat(paste0('    #enhancers after selecting N: ', length(enhancers), '\n'))
    }

    mcols(enhancers)$type <- as.factor('enhancer')
    enhancers
}
