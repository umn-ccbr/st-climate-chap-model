"""
prepare_chap_data.py
====================
Transforms data/harmonized_data_filled.csv into data/harmonized_data_chap.csv,
which uses the column names that CHAP (chap-core) expects:

  time_period  – YYYY-MM  (was: time = YYYYMM)
  location     – orgunit name (was: orgunitname)
  population   – population  (was: pop)
  preci        – precipitation (unchanged)
  temp_max     – maximum temperature (unchanged)
  disease_cases– disease cases, rounded to integer (CARBayesST requires int)

Every location is expanded to the global time range so CHAP sees no gaps
and never needs to NaN-pad covariates when it splits train/test windows.
disease_cases is left as NaN for genuinely missing months (CARBayesST
treats them as missing-Y and imputes during MCMC).

Usage:
    python3 scripts/prepare_chap_data.py \
        [--input data/harmonized_data_filled.csv] \
        [--output data/harmonized_data_chap.csv]
"""

import argparse
from pathlib import Path

import pandas as pd


def prepare(input_path: str, output_path: str) -> None:
    df = pd.read_csv(input_path)

    # ---------- column renames ----------
    rename_map: dict[str, str] = {}
    if "orgunitname" in df.columns:
        rename_map["orgunitname"] = "location"
    if "pop" in df.columns:
        rename_map["pop"] = "population"
    df = df.rename(columns=rename_map)

    # ---------- time conversion: YYYYMM → YYYY-MM ----------
    if "time" in df.columns and "time_period" not in df.columns:
        t = df["time"].astype(str).str.strip()
        df["time_period"] = t.str[:4] + "-" + t.str[4:6]
        df = df.drop(columns=["time"])

    # ---------- expand every location to the global monthly range ----------
    # CHAP aligns all locations to the global time range when it loads the
    # dataset, padding missing rows with NaN.  To prevent NaN covariates
    # reaching train.r, we pad here and ffill/bfill climate + population.
    # disease_cases is intentionally left NaN for missing months so that
    # CARBayesST treats them as missing-Y (imputed during MCMC).
    df["time_period"] = pd.PeriodIndex(df["time_period"], freq="M")
    df = df.sort_values(["location", "time_period"])

    global_start = df["time_period"].min()
    global_end   = df["time_period"].max()
    full_index   = pd.period_range(global_start, global_end, freq="M")

    covariate_cols = [c for c in ["preci", "temp_max", "population"] if c in df.columns]

    filled_parts: list[pd.DataFrame] = []
    for loc, g in df.groupby("location"):
        g = g.set_index("time_period").reindex(full_index)
        g["location"] = loc
        g.index.name = "time_period"
        # ffill then bfill covariates so padded rows are never NaN
        for col in covariate_cols:
            g[col] = g[col].ffill().bfill()
        filled_parts.append(g.reset_index())

    df = pd.concat(filled_parts, ignore_index=True)

    # Any location that was all-NaN for population gets the dataset median
    if "population" in df.columns:
        global_pop_median = df["population"].median()
        df["population"] = df["population"].fillna(global_pop_median)

    df["time_period"] = df["time_period"].astype(str)   # back to "YYYY-MM" string

    # ---------- round disease_cases to integer ----------
    # CARBayesST Poisson requires integer counts.  We use nullable Int64 so
    # genuinely missing months remain NaN rather than becoming 0.
    if "disease_cases" in df.columns:
        df["disease_cases"] = df["disease_cases"].round(0).astype("Int64")

    # ---------- column order ----------
    core_cols = ["time_period", "location", "population", "preci", "temp_max", "disease_cases"]
    extra_cols = [c for c in df.columns if c not in core_cols + ["orgunitid"]]
    out_cols   = [c for c in core_cols + extra_cols if c in df.columns]
    df = df[out_cols]

    Path(output_path).parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(output_path, index=False)

    n_pop_null = df["population"].isna().sum() if "population" in df.columns else 0
    n_dc_null  = df["disease_cases"].isna().sum() if "disease_cases" in df.columns else 0
    print(f"Written {len(df):,} rows → {output_path}")
    print(f"  Columns      : {list(df.columns)}")
    print(f"  time_period  : {df['time_period'].iloc[0]} … {df['time_period'].iloc[-1]}")
    print(f"  locations    : {df['location'].nunique()}")
    print(f"  pop nulls    : {n_pop_null}")
    print(f"  disease nulls: {n_dc_null}  (expected — missing Y for MCMC imputation)")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Transform harmonized CSV to CHAP format")
    parser.add_argument("--input",  default="data/harmonized_data_filled.csv")
    parser.add_argument("--output", default="data/harmonized_data_chap.csv")
    args = parser.parse_args()
    prepare(args.input, args.output)
