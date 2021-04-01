library(GenomicRanges)
library(GenomicAlignments)
library(SummarizedExperiment)

# Use the histone reads to create profiles for a list of genomic sites
# The sites could for example be enhancer sites found with `find_enhancers`
create_profiles <- function(ranges, bam_file, reference, name = NA, bin_size = 100, ignore_strand = FALSE)
{
    window_size <- unique(width(ranges))
    if (length(window_size) > 1) {
        stop('Not all ranges are of the same length.')
    }

    # Explode the ranges into bins.
    # We want each bin to be of a guarenteed size (bin_size). This means we
    # might have to skip some data when a range is not an exact multiple of
    # bin_size. However, the tile() function does not allow for this, so we
    # need a bit of a hack with slidingWindows.
    ranges_binned <- slidingWindows(ranges, width = bin_size, step = bin_size)
    ranges_binned <- ranges_binned[width(ranges_binned) == bin_size]  # Remove incomplete bins

    num_reads <- countBam(bam_file)$records

    # Compute the number of reads that overlap between the bam_file and the
    # binned ranges.
    if (is.na(Rsamtools::yieldSize(bam_file))) {
        # Read entire BAM file in one chunk
        histone <- unlist(readGAlignmentsList(bam_file))
        histone <- histone[from(findOverlaps(histone, ranges_binned, ignore.strand = ignore_strand))]
        profiles <- assays(
            summarizeOverlaps(
                unlist(ranges_binned),
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
            histone <- histone[from(findOverlaps(histone, ranges_binned, ignore.strand = ignore_strand))]
            counts <- assays(
                summarizeOverlaps(
                    unlist(ranges_binned),
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
    profiles <- matrix(profiles, nrow = length(ranges), byrow = TRUE)  # ranges x bins
    attr(profiles, 'ranges') <- ranges
    attr(profiles, 'window_size') <- window_size 
    attr(profiles, 'bin_size') <- bin_size 
    attr(profiles, 'num_reads') <- num_reads  # Needed when normalizing profiles

    if (is.na(name)) {
        name <- tools::file_path_sans_ext(basename(bam_file))
    }
    colnames(profiles) <- rep(name, ncol(profiles))

    if (!ignore_strand) {
        to_reverse <- as.vector(strand(ranges) == "-")
        profiles[to_reverse,] <- profiles[to_reverse, ncol(profiles):1]
    }

    if (!is.null(reference)) {
        profiles <- normalize_profiles(profiles, reference)
    }

    class(profiles) <- c('Profiles', 'matrix')
    profiles
}

# Normalize one set of profiles with another
normalize_profiles <- function(profiles, reference) {
    scaling_factor <- attr(profiles, 'num_reads') / attr(reference, 'num_reads')
    profiles <- profiles - scaling_factor * reference

    # Make proper counts (positive integer values)
    round(pmax(profiles, 0))
}


# Method functions for Profiles
###############################
is.Profiles <- function(x) inherits(x, 'Profiles')

# Get the type (enhancer, promoter or random) of each profile
profile_type <- function(x, ...) UseMethod("profile_type", x)
profile_type.Profiles <- function(x) attr(x, 'ranges')$type

# Interact with use the data part of the Profiles object
`profile_data<-` <- function(x, ...) UseMethod("profile_data<-", x)
`profile_data<-.Profiles` <- function(profiles, values)
{
    attr(values, 'ranges') <- attr(profiles, 'ranges')
    attr(values, 'window_size') <- attr(profiles, 'window_size')
    attr(values, 'bin_size') <- attr(profiles, 'bin_size')
    attr(values, 'num_reads') <- attr(profiles, 'num_reads')
    class(values) <- c('Profiles', 'matrix') 
    values
}

profile_data <- function(x, ...) UseMethod("profile_data<-", x)
profile_data.Profiles <- function(profiles)
{
    as.matrix(profiles)
}

# When subsetting a Profiles object, keep the metadata aligned
`[.Profiles` <- function(profile, i, ...) {
    selection <- NextMethod()
    attr(selection, 'num_reads') <- attr(profile, 'num_reads')
    attr(selection, 'window_size') <- attr(profile, 'window_size')
    attr(selection, 'bin_size') <- attr(profile, 'bin_size')
    if (missing(i)) {
        attr(selection, 'ranges') <- attr(profile, 'ranges')
    } else {
        attr(selection, 'ranges') <- attr(profile, 'ranges')[i]
    }
    class(selection) <- c('Profiles', 'matrix')
    selection
}

as.matrix.Profiles <- function(profiles)
{
    attr(profiles, 'ranges') <- NULL
    attr(profiles, 'window_size') <- NULL
    attr(profiles, 'bin_size') <- NULL
    attr(profiles, 'num_reads') <- NULL
    class(profiles) <- 'matrix'
    profiles
}

# Concatenate two or more Profiles row-wise
rbind.Profiles <- function(...)
{
    result <- do.call(rbind, lapply(list(...), as.matrix))
    attr(result, 'ranges') <- do.call(c, sapply(list(...), function(x) attr(x, 'ranges')))
    attr(result, 'num_reads') <- sum(sapply(list(...), function(x) attr(x, 'num_reads')))
    attr(result, 'window_size') <- attr(list(...)[[1]], 'window_size')
    attr(result, 'bin_size') <- attr(list(...)[[1]], 'bin_size')
    class(result) <- c('Profiles', 'matrix')
    result
}

# Concatenate two or more Profiles column-wise
cbind.Profiles <- function(...)
{
    result <- do.call(cbind, lapply(list(...), as.matrix))
    attr(result, 'ranges') <- attr(list(...)[[1]], 'ranges')
    attr(result, 'num_reads') <- attr(list(...)[[1]], 'num_reads')
    attr(result, 'window_size') <- attr(list(...)[[1]], 'window_size')
    attr(result, 'bin_size') <- attr(list(...)[[1]], 'bin_size')
    class(result) <- c('Profiles', 'matrix')
    result
}

names.Profiles <- function(profiles)
{
    unique(colnames(profiles))
}
