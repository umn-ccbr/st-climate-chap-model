#!/usr/bin/env python3
import sys
import pandas as pd

def main():
    if len(sys.argv) != 4:
        print("Usage: python split.py <input.csv> <train_out.csv> <test_out.csv>")
        sys.exit(1)

    input_path = sys.argv[1]
    train_out = sys.argv[2]
    test_out = sys.argv[3]

    # Load
    df = pd.read_csv(input_path)

    # Ensure time is integer-like (YYYYMM)
    df["time"] = pd.to_numeric(df["time"], errors="coerce")

    # Define boundaries
    TRAIN_MAX = 202412
    TEST_MIN = 202501

    train_df = df[(df["time"] >= 201701) & (df["time"] <= TRAIN_MAX)]
    test_df  = df[df["time"] >= TEST_MIN]

    # Save
    train_df.to_csv(train_out, index=False)
    test_df.to_csv(test_out, index=False)

    print(f"Train rows: {len(train_df)} → {train_out}")
    print(f"Test rows: {len(test_df)} → {test_out}")

if __name__ == "__main__":
    main()
