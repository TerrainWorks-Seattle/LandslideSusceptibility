task <- mlr3spatiotempcv::TaskClassifST$new(
  "landslide",
  backend = mlr3::as_data_backend(mlr3spatiotempcv::ecuador),
  target = "slides",
  positive = "TRUE",
  extra_args = list(coordinate_names = c("x", "y"), crs = 32717)
)

learner <- mlr3::lrn(
  "classif.rpart",
  maxdepth = 3,
  predict_type = "prob"
)

resampling <- mlr3::rsmp("repeated_spcv_coords", folds = 4, repeats = 2)

result <- mlr3::resample(
  task = task,
  learner = learner,
  resampling = resampling
)

result$resampling$test_set(9)

result$predictions()

result$aggregate(measures = mlr3::msr("classif.ce"))

mlr3spatiotempcv::autoplot(resamplingMethod, task, fold_id = c(1:2))
