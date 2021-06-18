library(yaml)
library(argparser)
library(caret)
library(preprint)

config <- read_yaml('workflow/config.yaml')

# parser <- arg_parser('Run the train-predict loop')
# parser <- add_argument(parser, 'cell_line', type = 'character',
#                        help = paste0('The cell line to process. Either ', paste(config$cell_lines, collapse = ' or ')))
# cell_line <- parse_args(parser)$cell_line
cell_line <- 'K562'

# Load the profiles
profiles <- readRDS(paste0(config$data_dir, '/', cell_line, '/data_R/profiles.rds'))

# Create cross validation groups
folds <- caret::createFolds(profile_type(profiles), k = 5)

# Evaluate performance using 5-fold cross validation, and parameter tuning
# using an inner 5-fold cross-validation loop.
predictions = list()
reference = list()
for (fold in folds) {
    characteristic_profiles <- aggregate_profiles(profiles[-fold,])
    train_data <- pattern_likelihoods(profiles[-fold,], characteristic_profiles, measure = 'Bayesian')
    test_data <- pattern_likelihoods(profiles[fold,], characteristic_profiles, measure = 'Bayesian')

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

fname <- paste0(config$data_dir, '/', cell_line, '/data_R/predictions.RData')
dir.create(dirname(fname), recursive = TRUE, showWarnings = FALSE)
save(predictions, reference, file = fname)
cat(paste0('Saved predictions to ', fname, '\n'))
