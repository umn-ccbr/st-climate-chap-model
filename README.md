# ST-CAR spatiotemporal forecasting model (R)

This repository contains an R implementation of a spatiotemporal disease forecasting model intended as a CHAP-compatible example for integration with DHIS2/CHAP workflows. The core model is a Poisson spatiotemporal CAR model (CARBayesST::ST.CARar) that models area-level monthly disease counts with spatial random effects and temporal autoregression.

Key points (short):
- Modeling approach: ST.CARar (CARBayesST) — Poisson outcome with area-level spatial random effects and AR temporal structure.
- Target: area-month disease counts (variable `disease_cases` / `marlaria` in the pipeline).
- Main covariates: population offset (log(pop)), precipitation (preci → PRCP) with lags (lag1, lag2), and maximum temperature (temp_max → TEMPmax) with lags (lag1, lag2, lag3). Natural splines are used for most lagged climate covariates.
- Input/outputs: training uses `data/harmonized_data.csv` and outputs fitted model RDS plus quartile summaries and training metadata in `output/`.

# Test run with sample Malaria dataset

Here's a quick command you can run to start this model in its
container, and run an evaluation with the `chap` CLI tool:

```sh
$ docker compose up --build
$ chap eval --model-name 'http://localhost:9090' --dataset-csv ./test-run/test-run.csv --output-file ./test-run/test-run-eval.nc --run-config.is-chapkit-model
```

# Adding as a CHAP model

You can include the following line in `config/configured_models/default.yaml`:

```yaml
- url: https://github.com/umn-ccbr/st-climate-chap-model
  versions:
    main: "@main"
```

*This is currently in testing*, as I haven't been able to get this
model to appear by doing this in my own CHAP-core instance.

## Files 
- `train.r` — main training script. Loads augmented data from `model_helpers.R`, selects an optimal orgunit/time configuration from `output/tradeoff_orgunits_vs_months.csv` (produced by `scripts/find_complete_timespan.py`), subsets the adjacency matrix `W_orgunits_CARBayesST.rds`, computes spline knot quartiles, fits `ST.CARar`, and saves the fitted model (`.rds`) and quartiles.
- `model_helpers.R` — helper functions:
  - `augment_harmonized_data()` reads `data/harmonized_data.csv`, creates `spaceid`/`timeid`, centers climate predictors, and creates lagged covariates: `lag1_PRCP`, `lag2_PRCP`, `lag1_TEMPmax`, `lag2_TEMPmax`, `lag3_TEMPmax`.
  - `compute_lag_quartiles()` computes 25/50/75% quartiles for each lag variable and returns knot locations for splines.
- `predict.r` — prediction/inference utilities. Contains a (work-in-progress) MCMC prediction routine that re-uses training components (summary, design matrix, adjacency) to produce predictive samples. Also includes `predict_chap()` — a minimal wrapper that reads a saved model and future climate CSV and writes predictions.
- `scripts/find_complete_timespan.py` — helper to find recommended orgunits and time spans (used by training).
- `W_orgunits_CARBayesST.rds` — adjacency matrix keyed by `orgunitname` used for CAR spatial structure.
- `TRAINING_REFACTOR_SUMMARY.md` — details of the refactored training pipeline and the chosen optimal configuration.

## Model covariates and offsets
- Outcome: `disease_cases` (also referenced as `marlaria` in some helper code).
- Offset: log(pop) where `pop` is population per area-month.
- Climate covariates (centered across dataset):
  - Precipitation (`preci` → `PRCP`), used with lags: `lag1_PRCP`, `lag2_PRCP` (natural splines).
  - Maximum temperature (`temp_max` → `TEMPmax`), used as `TEMPmax` and lagged spline terms `lag1_TEMPmax`, `lag2_TEMPmax`, `lag3_TEMPmax`. (Note: a spline on `TEMPmax` was commented out in `train.r` due to collinearity; a linear temp term is included.)

## How training works (high level)
1. Data augmentation: `augment_harmonized_data()` constructs lagged covariates and indexing fields required by CARBayesST.
2. Configuration selection: `train.r` loads a tradeoff table and recommended orgunits to pick a balanced panel (orgunits × months) that maximizes usable observations and ends at the most recent month.
3. Adjacency handling: `W` is subset to selected orgunits. Isolated areas (zero neighbors) are given 1–3 random neighbors as a stopgap to avoid modeling issues.
4. Spline knots: quartiles of lag covariates are computed to place natural spline knots.
5. Fit: `ST.CARar(...)` with Poisson family, offset(log(pop)), spline/transformed covariates, AR(2) temporal structure, and spatial CAR prior.

## Prediction / outputs
- Fitted model saved as an RDS file (path passed to `predict_chap()` or other prediction utilities).
- Quartile/knot information saved alongside the model (model_fn.quartiles) to ensure consistent pre-processing at prediction time.
- Prediction utilities in `predict.r` produce MCMC draws and pointwise mean predictions; they expect the same covariate names and indexing as in training.

## CHAP and DHIS2 integration
- CHAP integration: The repository follows the CHAP model layout and provides MLproject/entry points in the original example. CHAP runs the `train` and `predict` entry points inside a Docker image. To integrate with CHAP ensure the MLproject entry points point to `train.r` and `predict.r`, and provide a Docker image with R and required packages (`CARBayesST`, `splines`, and dependencies).
- DHIS2 integration / operationalization: The model expects area-level monthly inputs keyed by `orgunitname` and `time` and will output area-month forecasts. For DHIS2 integration you can:
  - Export DHIS2 aggregated case counts and population for each org unit and month into the harmonized CSV format used here.
  - Run the CHAP evaluation or the prediction wrapper to produce forecast CSVs, then map forecast rows back to DHIS2 orgunit IDs and time periods and push as aggregated data values or events via the DHIS2 API.
  - For production, automate harmonization (DHIS2→CSV) and prediction schedule, and implement a small adapter to convert CSV outputs into DHIS2-compatible payloads.
