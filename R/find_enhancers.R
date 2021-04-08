#' Find enhancer sites suitable to serve as training data.
#' @importFrom IRanges end findOverlaps from resize start start<- width width<-
#' @importFrom GenomicRanges strand<- seqnames seqinfo mcols mcols<-
#' @importFrom GenomeInfoDb seqlengths seqlengths<-
#' @export
find_enhancers  <- function(p300, DNase, window = 1000, N = NULL,
                            TSS = NULL, min_dist_to_promoter = 2000,
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

    if (!is.null(TSS)) { 
        # Remove p300 peaks that are too close to a promotor range, as specified by
        # min_dist_to_promoter.
        # We ignore the strand, which means all strands are presumed to be '+'.
        dist <- GenomicRanges::distanceToNearest(p300, TSS, ignore.strand = TRUE)
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
