View(mlr3spatiotempcv::ecuador)

task <- mlr3spatiotempcv::TaskClassifST$new(
  "landslide",
  backend = mlr3::as_data_backend(mlr3spatiotempcv::ecuador),
  target = "slides",
  positive = "TRUE",
  extra_args = list(coordinate_names = c("x", "y"), crs = 32717)
)

# Define resampling method for repeated k-fold spatial cross-validation
resampling <- mlr3::rsmp("repeated_spcv_coords", folds = 4, repeats = 10)

mlr3spatiotempcv::autoplot(resampling, task, fold_id = c(1:4), repeats_id = 1)

# Define a random forest learner
learner <- mlr3::lrn("classif.randomForest")

# Perform resampling method with the random forest learner and the landslides task
result <- mlr3::resample(
  task = task,
  learner = learner,
  resampling = resampling
)

# Display ROC plot for iteration 1
mlr3viz::autoplot(result$predictions()[[1]], type = "roc")

# Rows used for testing in iteration 1
result$resampling$test_set(1)

# Predictions for iteration 1
result$predictions()[[1]]

# Calculate average classification error
result$aggregate(measures = mlr3::msr("classif.ce"))

# Calculate average area under ROC curve
result$aggregate(measures = mlr3::msr("classif.auc"))

# Calculate average accuracy (proportion of correct predictions out of all predictions)
result$aggregate(measures = mlr3::msr("classif.acc"))
