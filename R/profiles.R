#' Create profiles for a list of genome ranges, using the genome coverage from a BAM file.
#'
#' @param ranges
#'     The ranges for which to make profiles. These should be windows of a
#'     fixed size, centered on a site of interest. These ranges can for example
#'     be created by the find_enhancers, find_promoters or find_random
#'     functions.
#' @param bam_file
#'     The BAM File to read the coverage data from.
#' @param reference
#'     The profiles to use as reference. The coverage counts from the reference
#'     are subtracted from the profile. Can be NULL to specify no reference.
#' @param name
#'     A name for the features which coverage we are computing. This is used to
#'     build the colnames() for the profiles. Defaults to NA, meaning the name
#'     will be automatically inferred from the BAM file.
#' @param bin_size
#'     The width of the bins, i.e. the resolution of the profiles.
#' @param ignore_strand
#'     Whether to ignore the strand when computing coverage.
#'
#' @return The profiles.
#'
#' @importFrom IRanges from width slidingWindows findOverlaps
#' @export
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

    num_reads <- Rsamtools::countBam(bam_file)$records

    # Compute the number of reads that overlap between the bam_file and the
    # binned ranges.
    if (is.na(Rsamtools::yieldSize(bam_file))) {
        # Read entire BAM file in one chunk
        histone <- unlist(GenomicAlignments::readGAlignmentsList(bam_file))
        histone <- histone[from(findOverlaps(histone, ranges_binned, ignore.strand = ignore_strand))]
        profiles <- SummarizedExperiment::assays(
            GenomicAlignments::summarizeOverlaps(
                unlist(ranges_binned),
                histone,
                inter.feature = FALSE,
                ignore.strand = ignore_strand
            )
        )$counts
    } else {
        # Read BAM file in chunks. Compute overlaps for each chunk.
        pb <- utils::txtProgressBar(min = 0, max = ceiling(num_reads / Rsamtools::yieldSize(bam_file)),
                                    style = 3)
        if (!Rsamtools::isOpen(bam_file)) {
            open(bam_file)
            on.exit(Rsamtools::close(bam_file))
            on.exit(close(pb))
        }
        profiles <- NULL
        utils::setTxtProgressBar(pb, 1)
        while (length(alignments <- GenomicAlignments::readGAlignmentsList(bam_file))) {
            histone = unlist(alignments)
            histone <- histone[from(findOverlaps(histone, ranges_binned, ignore.strand = ignore_strand))]
            counts <- SummarizedExperiment::assays(
                GenomicAlignments::summarizeOverlaps(
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
            utils::setTxtProgressBar(pb, utils::getTxtProgressBar(pb) + 1)
        }
    }

    profiles <- matrix(profiles, nrow = length(ranges), byrow = TRUE)  # ranges x bins

    if (is.na(name)) {
        name <- tools::file_path_sans_ext(basename(bam_file$path))
    }
    
    # name.1, name.2, name.3, name.4, ...
    colnames(profiles) <- lapply(seq(ncol(profiles)), function (x) paste0(name, '.', x))

    # Annotate with metadata
    attr(profiles, 'ranges') <- ranges
    attr(profiles, 'window_size') <- window_size 
    attr(profiles, 'bin_size') <- bin_size 
    attr(profiles, 'num_reads') <- rep(num_reads, ncol(profiles))  # Needed when normalizing profiles
    names(attr(profiles, 'num_reads')) <- colnames(profiles)

    if (!ignore_strand) {
        to_reverse <- as.vector(GenomicRanges::strand(ranges) == "-")
        profiles[to_reverse,] <- profiles[to_reverse, ncol(profiles):1]
    }

    if (!is.null(reference)) {
        profiles <- subtract_reference(profiles, reference)
    }

    class(profiles) <- c('Profiles', 'matrix')
    profiles
}

# Subtract one set of profiles from another
subtract_reference <- function(profiles, profiles_to_subtract) {
    scaling_factor <- attr(profiles, 'num_reads') / attr(profiles_to_subtract, 'num_reads')
    profiles <- profiles - scaling_factor * profiles_to_subtract

    # Make proper counts (positive integer values)
    profiles[profiles < 0] <- 0
    round(profiles)
}

#' Normalize one set of profiles with another
#' @export
normalize_profiles <- function(profiles, reference) {
    # Make sure the profiles are in the same order
    if (!identical(colnames(profiles), colnames(reference)))
        reference <- reference[, colnames(profiles)]

    scaling_factor <- attr(reference, 'num_reads') / attr(profiles, 'num_reads')
    round(profiles * scaling_factor)
}


# Method functions for Profiles
###############################

#' Check whether a given object is a Profiles object (or subclass).
#' @export
is.Profiles <- function(x) inherits(x, 'Profiles')

#' Get the type (enhancer, promoter or random) of each profile
#' @export
profile_type <- function(x, ...) UseMethod("profile_type", x)

#' Get the type (enhancer, promoter or random) of each profile
#' @export
profile_type.Profiles <- function(x) attr(x, 'ranges')$type

#' Interact with use the data part of the Profiles object
#' @export
`profile_data<-` <- function(x, ...) UseMethod('profile_data<-', x)

#' Interact with use the data part of the Profiles object
#' @export
`profile_data<-.Profiles` <- function(profiles, value)
{
    attr(value, 'ranges') <- attr(profiles, 'ranges')
    attr(value, 'window_size') <- attr(profiles, 'window_size')
    attr(value, 'bin_size') <- attr(profiles, 'bin_size')
    attr(value, 'num_reads') <- attr(profiles, 'num_reads')
    class(value) <- c('Profiles', 'matrix') 
    value
}

#' Get the data matrix of a Profiles object.
#' @export
profile_data <- function(x, ...) UseMethod('profile_data', x)

#' Get the data matrix of a Profiles object.
#' @export
profile_data.Profiles <- function(profiles)
{
    as.matrix(profiles)
}

#' When subsetting a Profiles object, keep the metadata aligned
#' @export
`[.Profiles` <- function(profile, i, j, ...) {
    selection <- NextMethod()
    attr(selection, 'window_size') <- attr(profile, 'window_size')
    attr(selection, 'bin_size') <- attr(profile, 'bin_size')
    if (missing(i)) {
        attr(selection, 'ranges') <- attr(profile, 'ranges')
    } else {
        attr(selection, 'ranges') <- attr(profile, 'ranges')[i]
    }
    if (missing(j)) {
        attr(selection, 'num_reads') <- attr(profile, 'num_reads')
    } else {
        attr(selection, 'num_reads') <- attr(profile, 'num_reads')[j]
    }
    class(selection) <- c('Profiles', 'matrix')
    selection
}

#' @export
as.matrix.Profiles <- function(profiles)
{
    attr(profiles, 'ranges') <- NULL
    attr(profiles, 'window_size') <- NULL
    attr(profiles, 'bin_size') <- NULL
    attr(profiles, 'num_reads') <- NULL
    class(profiles) <- 'matrix'
    profiles
}

#' Concatenate two or more Profiles row-wise
#' @export
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

#' Concatenate two or more Profiles column-wise
#' @export
cbind.Profiles <- function(...)
{
    result <- do.call(cbind, lapply(list(...), as.matrix))
    attr(result, 'ranges') <- attr(list(...)[[1]], 'ranges')
    attr(result, 'num_reads') <- do.call(c, lapply(list(...), attr, 'num_reads'))
    names(attr(result, 'num_reads')) <- colnames(result)
    attr(result, 'window_size') <- attr(list(...)[[1]], 'window_size')
    attr(result, 'bin_size') <- attr(list(...)[[1]], 'bin_size')
    class(result) <- c('Profiles', 'matrix')
    result
}

#' Get the names of the Profiles
#' @export
names.Profiles <- function(profiles)
{
    unique(colnames(profiles))
}
