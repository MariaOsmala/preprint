#' Compute for a set of profiles, the likelihood of them being an enhancer site.
#'
#' @description
#' Checks each profile against a list of characteristic profiles. For each
#' characteristic profile, a likelihood score is computed that reflects how
#' likely it is that the profile is of the type (enhancer, promoter, random,
#' etc.) represented by the characteristic profile. This score can be computed
#' either in Bayesian fashion (which is the main method advocated by PREPRINT)
#' or maximum-likelihood.
#'
#' @param profiles
#'     The profiles to compute the likelihoods for.
#' @param characteristic_profiles
#' 	   The characteristic profiles for each site type (enhancer, promoter, random, etc.)
#' @param measure
#' 	   The likelihood measure to use. Can be either "Bayesian" for Bayesian
#' 	   likelihood (the default), or "ML" for maximum likelihood.
#' @param positive_class
#' 	   When using Bayesian likelihood, you need to specify the site type
#' 	   ("enhancer", "promoter", "random", etc.) that is regarded as the
#' 	   positive class, i.e. the type you are looking to detect.
#' 	   Defaults to "enhancer".
#'
#' @return likelihoods
#'     A copy of the given Profiles object, but with the profile data replaced
#'     with likelihood scores. For each profile, likelihood scores are given
#'     for each of the given characteristic profiles. The score reflects the
#'     likelihood that the profile is of the type described by the
#'     characteristic profile.
#'
#' @export
pattern_likelihoods <- function(profiles, characteristic_profiles, measure = 'Bayesian', positive_class = 'enhancer')
{
    likelihoods = vector()
    for (name in names(profiles)) {
        sel <- colnames(profiles) == name

        if (measure == 'Bayesian') {
            pos_char_prof <- characteristic_profiles[[positive_class]]
            gamma <- estimate_gamma(profiles[, sel], pos_char_prof[sel])
            for (char_profile in characteristic_profiles) {
                likelihood <- likelihood_Bayesian(profiles[, sel], char_profile[sel], gamma)
                likelihoods <- cbind(likelihoods, likelihood)
            }
        } else if (measure == 'ML') {
            for (char_profile in characteristic_profiles) {
                likelihood <- likelihood_ML(profiles[, sel], char_profile[sel])
                likelihoods <- cbind(likelihoods, likelihood)
            }
        } else {
            stop(paste0('Invalid measure "', measure, '". Must be either "Bayesian" or "ML"'))
        }
    }

    colnames(likelihoods) <- levels(interaction(names(profiles), names(agg_patterns), lex.order = TRUE))
    profile_data(profiles) <- likelihoods
    profiles
}


likelihood_ML <- function(profiles, agg_pattern)
{
    # Formula 4 of the paper
    alpha <- rowSums(profiles) / sum(agg_pattern)

    # Formula 10 of the paper.
    lambda <- outer(alpha, agg_pattern)
    likelihood <- rowSums(stats::dpois(x = profiles, lambda = lambda, log = TRUE))
}

estimate_gamma <- function(profiles, pos_agg_pattern)
{
    # Formula 4 of the paper
    alpha <- rowSums(profiles) / sum(pos_agg_pattern)

    # Estimate gamma
    gamma <- try(MASS::fitdistr(alpha[alpha > 0], 'gamma', lower = list(shape = 1E-3, rate = 1E-3), upper = list(shape = 1E3, rate = 1E3), start = list(shape=1, rate=1))$estimate)
    if (class(gamma) == "try-error") {
        cat('Scaling down and retrying...')
        gamma <- MASS::fitdistr(alpha[alpha > 0] / 1000, 'gamma', lower = list(shape = 1E-3, rate = 1E-3), upper = list(shape = 1E3, rate = 1E3), start = list(shape=1, rate=1))$estimate
        gamma[['rate']] <- gamma[['rate']] / 1000
        cat(' success!\n')
    }

    gamma
}

likelihood_Bayesian <- function(profiles, agg_pattern, gamma)
{
    # Formula 5 of the paper.
    a0 <- gamma[['shape']]
    b0 <- gamma[['rate']]
    likelihood <- a0 * log(b0) + rowSums(profiles * log(agg_pattern)) - (a0 + rowSums(profiles)) * log(b0 + sum(agg_pattern)) + lgamma(a0 + rowSums(profiles)) - lgamma(a0) - rowSums(lgamma(profiles + 1))
}

#' Aggregate multiple profiles to create characteristic profiles.
#'
#' @description
#' Given a Profiles object, this function groups the profiles by site type
#' (enhancer, promoter, random, ...) and computes the mean profile for each type.
#' These aggregated profiles can then be used as characteristic profiles in the
#' pattern_likelihoods function.
#'
#' @param profiles
#'     A Profiles object containing multiple profiles, potentially belonging to
#'     different types.
#'
#' @return agg_profiles
#'     A list with for each profile type, the mean across the given profiles,
#'     i.e. the characteristic profile.
#'
#' @export
aggregate_profiles <- function(profiles)
{
    profile_types <- levels(profile_type(profiles))
    agg_profiles <- list()
    for (typ in profile_types) {
        agg_profiles[[typ]] <- colMeans(profiles[profile_type(profiles) == typ,])
    }
    agg_profiles
}
