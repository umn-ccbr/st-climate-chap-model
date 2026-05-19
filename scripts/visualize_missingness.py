#!/usr/bin/env python3
"""
Visualize covariate missingness for the recommended orgunit set with complete disease_cases.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

# Set style
plt.style.use('seaborn-v0_8-whitegrid')
plt.rcParams['figure.figsize'] = (14, 10)

# Load data
print("Loading data...")
df = pd.read_csv("data/harmonized_data_filled.csv")

# Drop current month
current_month = 202511
df = df[df['time'] != current_month].copy()

# Create timeid
sorted_times = sorted(df['time'].unique())
time_to_id = {t: i+1 for i, t in enumerate(sorted_times)}
df['timeid'] = df['time'].map(time_to_id)

# Filter for timeid > 4
df = df[df['timeid'] > 4].copy()

# Load the ORIGINAL recommended orgunits (based on disease_cases only)
# We want to visualize missingness in the 138 units with complete disease_cases
print("Loading original recommended orgunits (disease_cases complete)...")

# We need to find the original 138 units - let me recreate that analysis
# Get all units
all_units = sorted(df['orgunitname'].unique())
all_times = sorted(df['timeid'].unique())

# For each time, check disease_cases completeness
time_complete = []
for t in all_times:
    time_data = df[df['timeid'] == t]
    units_present = set(time_data['orgunitname'].unique())
    all_present = units_present == set(all_units)
    if all_present:
        cases_complete = time_data['disease_cases'].notna().all()
    else:
        cases_complete = False
    time_complete.append(all_present and cases_complete)

# Find units with best disease_cases coverage
unit_completeness = []
for unit in all_units:
    unit_data = df[df['orgunitname'] == unit]
    times_present = len(unit_data)
    cases_complete = unit_data['disease_cases'].notna().sum()
    unit_completeness.append({
        'orgunitname': unit,
        'times_present': times_present,
        'cases_complete': cases_complete
    })

unit_df = pd.DataFrame(unit_completeness)
sorted_units = unit_df.sort_values('cases_complete', ascending=False)

# Find the 138-unit configuration
def find_longest_recent_span(logical_array):
    if not logical_array.any():
        return {'start': None, 'end': None, 'length': 0}
    end_idx = len(logical_array) - 1
    length = 0
    for i in range(end_idx, -1, -1):
        if logical_array[i]:
            length += 1
        else:
            break
    if length > 0:
        return {'start': end_idx - length + 1, 'end': end_idx, 'length': length}
    true_indices = np.where(logical_array)[0]
    if len(true_indices) == 0:
        return {'start': None, 'end': None, 'length': 0}
    last_true = true_indices[-1]
    length = 1
    for i in range(last_true - 1, -1, -1):
        if logical_array[i]:
            length += 1
        else:
            break
    return {'start': last_true - length + 1, 'end': last_true, 'length': length}

# Find 138-unit config
for n_units in range(len(all_units), 100, -1):
    subset_units = set(sorted_units.head(n_units)['orgunitname'].tolist())
    subset_df = df[df['orgunitname'].isin(subset_units)]
    subset_complete = []
    for t in all_times:
        time_data = subset_df[subset_df['timeid'] == t]
        if len(time_data) != n_units:
            subset_complete.append(False)
            continue
        subset_complete.append(time_data['disease_cases'].notna().all())
    subset_complete = np.array(subset_complete)
    span = find_longest_recent_span(subset_complete)
    
    if span['length'] == 102 and all_times[span['end']] == 106:  # Ends at timeid 106
        selected_units = list(subset_units)
        start_timeid = all_times[span['start']]
        end_timeid = all_times[span['end']]
        print(f"Found configuration: {n_units} units, {span['length']} months")
        print(f"  Time range: {start_timeid} to {end_timeid}")
        break

# Filter to selected units and timespan
df_filtered = df[
    df['orgunitname'].isin(selected_units) &
    (df['timeid'] >= start_timeid) &
    (df['timeid'] <= end_timeid)
].copy()

print(f"\nFiltered data shape: {len(df_filtered)} rows")
print(f"Units: {len(selected_units)}")
print(f"Times: {len(df_filtered['timeid'].unique())}")

# Compute lagged variables
def add_lags(data, col, max_lag):
    for lag in range(1, max_lag + 1):
        lag_col = f'lag{lag}_{col}'
        data[lag_col] = np.nan
        for unit in data['orgunitname'].unique():
            mask = data['orgunitname'] == unit
            unit_data = data.loc[mask].sort_values('timeid')
            unit_indices = unit_data.index
            values = unit_data[col].values
            if len(values) > lag:
                data.loc[unit_indices[lag:], lag_col] = values[:-lag]
    return data

# Center climate variables
df_filtered['PRCP'] = df_filtered['preci'] - df_filtered['preci'].mean()
df_filtered['TEMPmax'] = df_filtered['temp_max'] - df_filtered['temp_max'].mean()

# Add lags
df_filtered = add_lags(df_filtered, 'PRCP', 2)
df_filtered = add_lags(df_filtered, 'TEMPmax', 3)

# Define model variables
model_vars = ['disease_cases', 'pop', 'preci', 'temp_max']

# Calculate missingness
print("\nCalculating missingness percentages...")
missingness = {}
for var in model_vars:
    missing_count = df_filtered[var].isna().sum()
    total_count = len(df_filtered)
    pct = (missing_count / total_count) * 100
    missingness[var] = pct
    print(f"  {var}: {pct:.1f}% missing ({missing_count}/{total_count})")

# Create visualization
fig, axes = plt.subplots(2, 2, figsize=(16, 12))
fig.suptitle(f'Covariate Missingness Analysis\n{len(selected_units)} Units × 102 Months (May 2017 - Oct 2025)\nComplete disease_cases, Sparse Covariates', 
             fontsize=16, fontweight='bold')

# 1. Bar chart of missingness percentages
ax1 = axes[0, 0]
vars_sorted = sorted(missingness.items(), key=lambda x: x[1], reverse=True)
var_names = [v[0] for v in vars_sorted]
var_pcts = [v[1] for v in vars_sorted]
colors = ['green' if p == 0 else 'orange' if p < 50 else 'red' for p in var_pcts]

bars = ax1.barh(var_names, var_pcts, color=colors, alpha=0.7, edgecolor='black')
ax1.set_xlabel('Percentage Missing (%)', fontsize=12)
ax1.set_title('Missingness by Variable', fontsize=14, fontweight='bold')
ax1.set_xlim(0, 100)
ax1.grid(axis='x', alpha=0.3)

# Add percentage labels
for i, (var, pct) in enumerate(vars_sorted):
    ax1.text(pct + 2, i, f'{pct:.1f}%', va='center', fontsize=10)

# 2. Missingness over time
ax2 = axes[0, 1]
time_missingness = []
for t in sorted(df_filtered['timeid'].unique()):
    time_data = df_filtered[df_filtered['timeid'] == t]
    pct = (time_data[model_vars].isna().sum().sum() / (len(time_data) * len(model_vars))) * 100
    time_missingness.append({'timeid': t, 'pct_missing': pct})

time_miss_df = pd.DataFrame(time_missingness)
ax2.plot(time_miss_df['timeid'], time_miss_df['pct_missing'], linewidth=2, color='darkred')
ax2.fill_between(time_miss_df['timeid'], time_miss_df['pct_missing'], alpha=0.3, color='red')
ax2.set_xlabel('Time ID', fontsize=12)
ax2.set_ylabel('Percentage Missing (%)', fontsize=12)
ax2.set_title('Overall Missingness Over Time', fontsize=14, fontweight='bold')
ax2.grid(alpha=0.3)
ax2.set_ylim(0, 100)

# 3. Heatmap of missingness by variable and time (sample)
ax3 = axes[1, 0]
# Sample every 5th time point for readability
sample_times = sorted(df_filtered['timeid'].unique())[::5]
heatmap_data = []
for t in sample_times:
    time_data = df_filtered[df_filtered['timeid'] == t]
    row = [time_data[var].isna().sum() / len(time_data) * 100 for var in model_vars]
    heatmap_data.append(row)

heatmap_df = pd.DataFrame(heatmap_data, columns=model_vars, index=sample_times)
im = ax3.imshow(heatmap_df.T, cmap='RdYlGn_r', aspect='auto', vmin=0, vmax=100)
ax3.set_xticks(range(len(sample_times)))
ax3.set_xticklabels(sample_times, rotation=45)
ax3.set_yticks(range(len(model_vars)))
ax3.set_yticklabels(model_vars)
cbar = plt.colorbar(im, ax=ax3)
cbar.set_label('% Missing', rotation=270, labelpad=20)
ax3.set_xlabel('Time ID (sampled every 5th)', fontsize=12)
ax3.set_ylabel('Variable', fontsize=12)
ax3.set_title('Missingness Heatmap (Time × Variable)', fontsize=14, fontweight='bold')

# 4. Summary statistics table
ax4 = axes[1, 1]
ax4.axis('off')

summary_text = f"""
SUMMARY STATISTICS
{'='*50}

Configuration:
  • Units: {len(selected_units)} (with complete disease_cases)
  • Time span: 102 months (May 2017 - Oct 2025)
  • Total observations: {len(df_filtered):,}

Missingness Overview:
  • Response (disease_cases): {missingness['disease_cases']:.1f}%
  • Population: {missingness['pop']:.1f}%
  • Climate variables: {np.mean([missingness['preci'], missingness['temp_max']]):.1f}% avg
  • Lagged variables: {np.mean([missingness[v] for v in model_vars if 'lag' in v]):.1f}% avg

Complete Cases:
  • Rows with all variables: {df_filtered[model_vars].notna().all(axis=1).sum():,}
  • Percentage complete: {(df_filtered[model_vars].notna().all(axis=1).sum() / len(df_filtered) * 100):.1f}%

Impact:
  • Original expected obs: {len(selected_units)} × 102 = {len(selected_units) * 102:,}
  • Actual obs: {len(df_filtered):,}
  • Complete obs (all vars): {df_filtered[model_vars].notna().all(axis=1).sum():,}
  • Loss: {len(selected_units) * 102 - df_filtered[model_vars].notna().all(axis=1).sum():,} obs ({((len(selected_units) * 102 - df_filtered[model_vars].notna().all(axis=1).sum()) / (len(selected_units) * 102) * 100):.1f}%)
"""

ax4.text(0.1, 0.9, summary_text, transform=ax4.transAxes, 
         fontsize=11, verticalalignment='top', fontfamily='monospace',
         bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.3))

plt.tight_layout()
output_path = Path("output/covariate_missingness.png")
output_path.parent.mkdir(exist_ok=True)
plt.savefig(output_path, dpi=300, bbox_inches='tight')
print(f"\n✓ Saved visualization to: {output_path}")

# Also create a simple CSV summary
summary_df = pd.DataFrame([
    {'variable': var, 'pct_missing': pct, 'n_missing': int(len(df_filtered) * pct / 100)}
    for var, pct in missingness.items()
])
summary_df = summary_df.sort_values('pct_missing', ascending=False)
summary_path = Path("output/covariate_missingness.csv")
summary_df.to_csv(summary_path, index=False)
print(f"✓ Saved summary to: {summary_path}")

plt.show()
