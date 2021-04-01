# Meat of the PREPRINT method
pattern_likelihoods <- function(profiles, agg_patterns, measure = 'ML', positive_class = 'enhancer')
{
    likelihoods = vector()
    for (name in names(profiles)) {
        sel <- colnames(profiles) == name

        if (measure == 'Bayesian') {
            pos_agg_pattern <- agg_patterns[[positive_class]]
            gamma <- estimate_gamma(profiles[, sel], pos_agg_pattern[sel])
            for (agg_pattern in agg_patterns) {
                likelihood <- likelihood_Bayesian(profiles[, sel], agg_pattern[sel], gamma)
                likelihoods <- cbind(likelihoods, likelihood)
            }
        } else if (measure == 'ML') {
            for (agg_pattern in agg_patterns) {
                likelihood <- likelihood_ML(profiles[, sel], agg_pattern[sel])
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
    likelihood <- rowSums(dpois(x = profiles, lambda = lambda, log = TRUE))
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

aggregate_patterns <- function(profiles)
{
    profile_types <- levels(profile_type(profiles))
    agg_patterns <- list()
    for (typ in profile_types) {
        agg_patterns[[typ]] <- colMeans(profiles[profile_type(profiles) == typ,])
    }
    agg_patterns
}
