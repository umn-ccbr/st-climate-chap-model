source("train.r")
source("predict.r")

train_chap(
  csv_fn      = "test-run/test-run.csv",
  model_fn    = "model.rds",
  polygons_fn = "test-run/test-run.geojson"
)

predict_chap(
  model_fn         = "model.rds",
  historic_data_fn = "test-run/test-run.csv",
  future_data_fn   = "test-run/test-run.csv",
  predictions_fn   = "test-run/predictions.csv"
)
