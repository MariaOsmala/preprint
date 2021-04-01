library(BSgenome.Hsapiens.UCSC.hg19)
library(GenomicRanges)

find_random <- function(window = 1000, N = NULL,
                        chroms_of_interest = NULL,
                        p300 = NULL, max_dist_to_p300 = 2000,
                        TSS = NULL, max_dist_to_TSS = 2000,
                        blacklist = NULL, verbose = TRUE)
{
    if (verbose) cat('Finding suitable random sites:\n')

    if (is.null(chroms_of_interest)) {
        chroms_of_interest = c('chr1',  'chr2',  'chr3',  'chr4',  'chr5',  'chr6',
                               'chr7',  'chr8',  'chr9', 'chr10', 'chr11', 'chr12',
                               'chr13', 'chr14', 'chr15', 'chr16', 'chr17',
                               'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX')
    }
    whole_genome <- GRanges(chroms_of_interest, IRanges(start=1, end=seqlengths(Hsapiens)[chroms_of_interest]))

    if (is.null(blacklist)) {
        blacklist = GRanges()
    }

    # Remove regions that are too close to a p300 peak, as specified by
    # max_dist_to_p300
    if (!is.null(p300)) {
        p300_exclusion_zone <- p300
        p300_exclusion_zone <- resize(p300_exclusion_zone, width = max_dist_to_p300 + (max_dist_to_p300 + 1) %% 2, fix='center')
        blacklist <- union(blacklist, p300_exclusion_zone)
        if (verbose) cat(paste0('    There are ', length(p300_exclusion_zone), ' p300 peaks to be avoided.\n'))
    }

    # Remove regions that are too close to an TSS annotation, as specified by
    # max_dist_to_TSS
    if (!is.null(TSS)) {
        TSS_exclusion_zone <- TSS
        TSS_exclusion_zone <- resize(TSS_exclusion_zone, width = max_dist_to_TSS + (max_dist_to_TSS + 1) %% 2, fix='center')
        blacklist <- union(blacklist, TSS_exclusion_zone)
        if (verbose) cat(paste0('    There are ', length(TSS_exclusion_zone), ' TSS annotated sites to be avoided.\n'))
    }

    blacklist <- keepSeqlevels(blacklist, chroms_of_interest, pruning.mode = 'coarse')

    # reduce() will greatly reduce the number of ranges by merging adjacent ranges.
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
