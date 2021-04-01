library(caret)

distance_measure="Bayes_estimated_priors"
path='~/scratch_cs/csb/projects/enhancer_prediction/aaltorse/Data'
cell_line='K562'

source('code/create_profiles.R')

# Load the profiles
load(paste0(path, "/", cell_line, "/data_R/profiles.RData"))

# Create cross validation groups
folds <- caret::createFolds(profile_type(profiles), k = k)

# Meat of the machine learning function
distance_to_summaries <- function(profiles, summaries, measure = 'ML')
{
    distances = vector()
    for (name in names(profiles)) {
        sel <- colnames(profiles) == name
        for (summary in summaries) {
            alpha <- rowSums(profiles[, sel]) / sum(summary[sel])
            lambdas <- outer(alpha, summary[sel])
            likelyhood <- rowSums(dpois(x = profiles[, sel], lambda = lambdas, log = TRUE))
            distances <- cbind(distances, likelyhood)
        }
    }
    colnames(distances) <- levels(interaction(names(profiles), names(summaries), lex.order = TRUE))
    profile_data(profiles) <- distances
    profiles
}

make_summaries <- function(profiles)
{
    profile_types <- levels(profile_type(profiles))
    summaries <- list()
    for (typ in profile_types) {
        summaries[[typ]] <- colMeans(profiles[profile_type(profiles) == typ,])
    }
    summaries
}

# Evaluate performance using 5-fold cross validation, and parameter tuning
# using an inner 5-fold cross-validation loop.
predictions = list()
reference = list()
for (fold in folds) {
    train_summaries <- make_summaries(profiles[-fold,])
    train_data <- distance_to_summaries(profiles[-fold,], train_summaries)
    test_data <- distance_to_summaries(profiles[fold,], train_summaries)

    train_labels <- as.factor(ifelse(profile_type(train_data) == "enhancer",
                                     "enhancer", "not.enhancer"))
    test_labels <- as.factor(ifelse(profile_type(test_data) == "enhancer",
                                    "enhancer", "not.enhancer"))

    tuneGrid <- expand.grid(C = 10 ^ seq(0, 2, length = 20), sigma = 1 / 45, Weight = 0.25)
    fitControl <- trainControl(method = "cv", number = 5, classProbs = TRUE, verboseIter = TRUE)
    model <- train(train_data, train_labels, method = 'svmRadialWeights',
                   trControl = fitControl, verbose = TRUE, tuneGrid = tuneGrid)
    pred <- predict(model, test_data)
    cat(paste0('Accuracy: ', sum(pred == test_labels) / length(test_labels), '\n'))

    predictions <- unlist(list(predictions, pred))
    reference <- unlist(list(reference, test_labels))
}
print(confusionMatrix(data = unlist(predictions), reference = unlist(reference)))
