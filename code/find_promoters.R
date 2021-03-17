library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg19)

find_promoters <- function(TSS_annotation, DNase, window = 1000, N = NULL,
                           between_TSS_distance = 2000, blacklist = NULL)
{
    # Find the closest overlap between DNase and TSS_annotation
    dist <- distanceToNearest(DNase, TSS_annotation, ignore.strand = TRUE)
    TSS_with_DNase <- TSS_annotation[to(dist)]
    mcols(TSS_with_DNase) <- mcols(DNase[from(dist)])[c("signalValue", "pValue")]
    mcols(TSS_with_DNase)$distance <- mcols(dist)$distance

    # Compute the distance from the peak in the DNase to the corresponding TSS_annotation
    DNase_peak <- DNase
    start(DNase_peak) <- start(DNase_peak) + DNase_peak$peak
    width(DNase_peak) <- 1
    dist_peak <- distanceToNearest(DNase_peak, TSS_annotation, ignore.strand = TRUE)
    mcols(TSS_with_DNase)$peakDistance <- mcols(dist_peak[from(dist)])$distance

    # Order
    TSS_with_DNase <- TSS_with_DNase[order(mcols(TSS_with_DNase)$distance, decreasing = FALSE)]

    # Promoter sites are TSS that overlap with DNase
    to_keep <- from(findOverlaps(TSS_with_DNase, DNase, ignore.strand = TRUE))
    promoters <- TSS_with_DNase[to_keep] #60528

    # Remove any TSS that have other TSS nearby
    promoters <- promoters[order(promoters$distance)]
    tmp <- resize(promoters, between_TSS_distance / 2, fix = "center")
    tmp <- findOverlaps(tmp, ignore.strand = TRUE)
    tmp <- tmp[from(tmp) != to(tmp)]  # Ranges always overlap themselves
    to_drop <- unique(c(from(tmp), to(tmp)))  # All ranges that overlap with any other ranges
    if (length(to_drop) > 0) {
        promoters <- promoters[-to_drop]
    }
    print(paste0("#promoters that are nicely isolated: ", length(promoters)))

    # ???
    promoters <- promoters[promoters$distance == 0]

    # Create a window around the promoters
    # If window%%2==0, add 1 to window, extend the enhancer window/2 downstream and window/2 upsteam
    # If window%%2==1, extend the enhancer window (window-1)/2 downstream and upstream
    promoters <- resize(promoters, width = window + (window + 1) %% 2 , fix = "center" )

    if (!is.null(blacklist)) {
        # Remove promoters that fall within blacklisted regions
        to_remove <- unique(from(findOverlaps(promoters, blacklist)))
        if (length(to_remove) > 0) {
            promoters <- promoters[-to_remove]
        }
        print(paste0("#promoters after blacklist: ", length(promoters)))
    }

    if (!is.null(N)) {
        # Sort first by qValue, then by signalValue, then by peakDistance
        promoters <- promoters[order(promoters$pValue, promoters$signalValue, decreasing = TRUE)]
        promoters <- promoters[order(promoters$peakDistance, decreasing = FALSE)]
        promoters <- promoters[1:N]
        print(paste0("#promoters after selecting N: ", length(promoters)))
    }

    promoters
}
