#' Find random sites in the genome that are not enhancers or promoters.
#'
#' @description
#' Find random sites in the genome that are suitable as training data. These
#' sites are guarenteed not to be potential enhancer or promoter sites.
#'
#' @param window
#'     The length of the window to cut around each site. This window can later
#'     be used to create the corresponding coverage profile using
#'     `create_profiles`. The window will be centered on the site.
#' @param N
#'     The number of promoter sites to find. Promoter sites with the highest
#'     qValue will be used (using signalValue to resolve ties).
#' @param chroms_of_interest
#'     Set this to a list of chromosome names to restrict the selection of
#'     random sites to these chromosomes.
#' @param p300
#'     A `GRanges` object containing the p300 peaks that will be used to
#'     indicate enhancer sites that will be avoided.
#' @param min_dist_to_enhancer
#'     Desired minimum distance to any enhancer sites.
#' @param TSS_annotations
#'     A `GRanges` object containing the TSS annotations that will be used to
#'     indicate promoter sites that will be avoided.
#' @param min_dist_to_promoter
#'     Desired minimum distance to any promoter sites.
#' @param blacklist
#'     An optional `GRanges` object containing areas that should be avoided
#'     when searching for suitable sites.
#' @param verbose
#'     A logical scalar indicating whether to print out informational messages
#'     during the search for suitable sites.
#'
#' @return
#'     A `GRanges` object containing windows centered around suitable random
#'     sites. The site itself lies in the center of the window. The `RGanges`
#'     object has a `type` metadata column that is set to `'random'`.
#'
#' @seealso [find_enhancers()] for finding suitable enhancer sites and
#'          [find_promoters()] for finding suitable promoter sites.
#' @export
#'
#' @importFrom IRanges IRanges resize reduce
#' @importFrom GenomicRanges mcols

find_random <- function(window = 1000, N = NULL,
                        chroms_of_interest = NULL,
                        p300 = NULL, max_dist_to_enhancer = 2000,
                        TSS_annotations = NULL, max_dist_to_promoter = 2000,
                        blacklist = NULL, verbose = TRUE)
{
    if (verbose) cat('Finding suitable random sites:\n')

    if (is.null(chroms_of_interest)) {
        chroms_of_interest = c('chr1',  'chr2',  'chr3',  'chr4',  'chr5',  'chr6',
                               'chr7',  'chr8',  'chr9', 'chr10', 'chr11', 'chr12',
                               'chr13', 'chr14', 'chr15', 'chr16', 'chr17',
                               'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX')
    }
    whole_genome <- GenomicRanges::GRanges(chroms_of_interest, IRanges(start=1, end=seqlengths(BSgenome.Hsapiens.UCSC.hg19::Hsapiens)[chroms_of_interest]))

    if (is.null(blacklist)) {
        blacklist = GenomicRanges::GRanges()
    }

    # Remove regions that are too close to a p300 peak, as specified by
    # max_dist_to_enhancer
    if (!is.null(p300)) {
        enhancer_exclusion_zone <- p300
        enhancer_exclusion_zone <- resize(enhancer_exclusion_zone, width = max_dist_to_enhancer + (max_dist_to_enhancer + 1) %% 2, fix='center')
        blacklist <- GenomicRanges::union(blacklist, enhancer_exclusion_zone)
        if (verbose) cat(paste0('    There are ', length(enhancer_exclusion_zone), ' enhancer sites to be avoided.\n'))
    }

    # Remove regions that are too close to a promoter, as specified by
    # max_dist_to_promoters
    if (!is.null(TSS_annotations)) {
        promoter_exclusion_zone <- TSS_annotations
        promoter_exclusion_zone <- resize(promoter_exclusion_zone, width = max_dist_to_promoter + (max_dist_to_promoter + 1) %% 2, fix='center')
        blacklist <- GenomicRanges::union(blacklist, promoter_exclusion_zone)
        if (verbose) cat(paste0('    There are ', length(promoter_exclusion_zone), ' promoter sites to be avoided.\n'))
    }

    blacklist <- GenomeInfoDb::keepSeqlevels(blacklist, chroms_of_interest, pruning.mode = 'coarse')

    # reduce() will greatly reduce the number of ranges by merging adjacent ranges.
    # This is important, as createRandomRegions will loop over them.
    blacklist <- reduce(blacklist)

    if (verbose) {
        cat(paste0('    Blacklisted ', length(blacklist), ' sites in total.\n'))
        cat('    Now finding random regions that satisfy the criteria...')
    }
    random_regions <- regioneR::createRandomRegions(N, length.mean = window, length.sd = 0, genome = whole_genome, mask = blacklist, non.overlapping = FALSE)
    random_regions <- sort(random_regions)
    mcols(random_regions)$type <- as.factor('random')

    if (verbose) cat(' done.\n')

    random_regions
}
