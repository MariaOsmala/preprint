library(GenomicRanges)
library(GenomicAlignments)
library(SummarizedExperiment)

# Use the histone reads to create a profile for a list of sites.
# The sites could for example be enhancer sites found with `find_enhancers`
create_profile <- function(sites, histone, reference, bin_size = 100)
{
    # Explode the sites into bins.
    # We want each bin to be of a guarenteed size (bin_size). This means we
    # might have to skip some data when a range is not an exact multiple of
    # bin_size. However, the tile() function does not allow for this, so we
    # need a bit of a hack with slidingWindows.
    sites <- slidingWindows(sites, width = bin_size, step = bin_size)
    sites <- sites[width(sites) == bin_size]  # Remove incomplete windows

    num_reads <- length(histone)

    # Compute the number of reads that overlap between the histone and the
    # binned sites.
    histone <- histone[from(findOverlaps(histone, sites))]
    data <- assays(summarizeOverlaps(unlist(sites), histone, inter.feature = FALSE))$counts
    data <- matrix(data, ncol = length(sites))  # nbins x nranges

    profile <- list(data = data, num_reads = num_reads)

    if (!is.null(reference)) {
        profile <- normalize_profile(profile, reference)
    }

    profile
}

# Normalize one profile with another one
normalize_profile <- function(profile, reference) {
    scaling_factor <- profile$num_reads / reference$num_reads
    profile$data <- profile$data - scaling_factor * reference$data
    profile
}
