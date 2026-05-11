source("train.r")
source("predict.r")

train_chap(
  "input/training_data_harmonized.csv",
  "output/model.bin"
)

predict_chap(
  "output/model.bin",
  "data/training_data_harmonized_filled.csv",
  "input/test_data_harmonized.csv",
  "output/predictions.csv"
)
