library(GenomicRanges)
library(GenomicAlignments)
library(SummarizedExperiment)

# Use the histone reads to create profiles for a list of sites.
# The sites could for example be enhancer sites found with `find_enhancers`
create_profiles <- function(sites, bam_file, reference, bin_size = 100, ignore_strand = FALSE)
{
    # Explode the sites into bins.
    # We want each bin to be of a guarenteed size (bin_size). This means we
    # might have to skip some data when a range is not an exact multiple of
    # bin_size. However, the tile() function does not allow for this, so we
    # need a bit of a hack with slidingWindows.
    sites_binned <- slidingWindows(sites, width = bin_size, step = bin_size)
    sites_binned <- sites_binned[width(sites_binned) == bin_size]  # Remove incomplete windows

    num_reads <- countBam(bam_file)$records

    # Compute the number of reads that overlap between the bam_file and the
    # binned sites.
    if (is.na(Rsamtools::yieldSize(bam_file))) {
        # Read entire BAM file in one chunk
        histone <- unlist(readGAlignmentsList(bam_file))
        histone <- histone[from(findOverlaps(histone, sites_binned, ignore.strand = ignore_strand))]
        profiles <- assays(
            summarizeOverlaps(
                unlist(sites_binned),
                histone,
                inter.feature = FALSE,
                ignore.strand = ignore_strand
            )
        )$counts
    } else {
        # Read BAM file in chunks. Compute overlaps for each chunk.
        pb <- txtProgressBar(min = 0, max = ceiling(num_reads / yieldSize(bam_file)),
                             style = 3)
        if (!isOpen(bam_file)) {
            open(bam_file)
            on.exit(close(bam_file))
            on.exit(close(pb))
        }
        profiles <- NULL
        setTxtProgressBar(pb, 0)
        while (length(alignments <- readGAlignmentsList(bam_file))) {
            histone = unlist(alignments)
            histone <- histone[from(findOverlaps(histone, sites_binned, ignore.strand = ignore_strand))]
            counts <- assays(
                summarizeOverlaps(
                    unlist(sites_binned),
                    histone,
                    inter.feature = FALSE,
                    ignore.strand = ignore_strand
                )
            )$counts
            if (is.null(profiles)) {
                profiles <- counts
            } else {
                profiles <- profiles + counts
            }
            setTxtProgressBar(pb, getTxtProgressBar(pb) + 1)
        }
    }

    # Assemble the profiles matrix, annotate with metadata
    profiles <- matrix(profiles, nrow = length(sites), byrow = TRUE)  # ranges x bins
    attr(profiles, 'ranges') <- sites_binned
    attr(profiles, 'num_reads') <- num_reads  # Needed when normalizing profiles

    if (!ignore_strand) {
        to_reverse <- as.vector(strand(sites) == "-")
        profiles[to_reverse,] <- profiles[to_reverse, ncol(profiles):1]
    }

    if (!is.null(reference)) {
        profiles <- normalize_profiles(profiles, reference)
    }

    profiles
}

# Normalize one set of profiles with another
normalize_profiles <- function(profiles, reference) {
    scaling_factor <- attr(profiles, 'num_reads') / attr(reference, 'num_reads')
    profiles <- profiles - scaling_factor * reference

    # Make proper counts (positive integer values)
    round(pmax(profiles, 0))
}
