library(caret)

path='~/scratch_cs/csb/projects/enhancer_prediction/aaltorse/Data'
cell_line='K562'

source('code/create_profiles.R')
source('code/statistics.R')

# Load the profiles
load(paste0(path, "/", cell_line, "/data_R/profiles.RData"))

# Create cross validation groups
folds <- caret::createFolds(profile_type(profiles), k = 5)

# Evaluate performance using 5-fold cross validation, and parameter tuning
# using an inner 5-fold cross-validation loop.
predictions = list()
reference = list()
for (fold in folds) {
    train_patterns <- aggregate_patterns(profiles[-fold,])
    train_data <- pattern_likelihoods(profiles[-fold,], train_patterns, measure = 'Bayesian')
    test_data <- pattern_likelihoods(profiles[fold,], train_patterns, measure = 'Bayesian')

    train_labels <- as.factor(ifelse(profile_type(train_data) == "enhancer",
                                     "enhancer", "not.enhancer"))
    test_labels <- as.factor(ifelse(profile_type(test_data) == "enhancer",
                                    "enhancer", "not.enhancer"))

    tuneGrid <- expand.grid(C = 2 ^ seq(-5, 8, length = 20), sigma = 1 / 45, Weight = 0.25)
    fitControl <- trainControl(method = "cv", number = 5, classProbs = TRUE, verboseIter = TRUE)
    model <- train(train_data, train_labels, method = 'svmRadialWeights',
                   trControl = fitControl, tuneGrid = tuneGrid)
    pred <- predict(model, test_data)
    print(model)
    cat(paste0('Accuracy: ', sum(pred == test_labels) / length(test_labels), '\n'))

    predictions <- unlist(list(predictions, pred))
    reference <- unlist(list(reference, test_labels))
}
print(confusionMatrix(data = unlist(predictions), reference = unlist(reference)))
