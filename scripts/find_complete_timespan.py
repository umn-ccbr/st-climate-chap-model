#!/usr/bin/env python3
"""
Find the longest time span where all spatial units have complete disease_cases data.
Analyzes data/harmonized_data.csv and reports the time span.
"""

import pandas as pd
import numpy as np
from pathlib import Path


def find_longest_consecutive_span(logical_array):
    """
    Find the longest consecutive sequence of True values in a boolean array.
    
    Returns:
        dict with 'start' (index), 'end' (index), and 'length' of longest span
    """
    if not logical_array.any():
        return {'start': None, 'end': None, 'length': 0}
    
    # Find runs of consecutive True values
    # Pad with False at boundaries to detect starts and ends
    padded = np.concatenate([[False], logical_array, [False]])
    diff = np.diff(padded.astype(int))
    
    starts = np.where(diff == 1)[0]
    ends = np.where(diff == -1)[0]
    
    if len(starts) == 0:
        return {'start': None, 'end': None, 'length': 0}
    
    lengths = ends - starts
    max_idx = np.argmax(lengths)
    
    return {
        'start': starts[max_idx],
        'end': ends[max_idx] - 1,  # inclusive end
        'length': lengths[max_idx]
    }


def main():
    # Hardcoded path
    data_path = Path("data/harmonized_data.csv")
    
    if not data_path.exists():
        print(f"Error: Data file not found at {data_path}")
        return
    
    print(f"Loading data from: {data_path}")
    df = pd.read_csv(data_path)
    
    # Check required columns - now including ALL variables needed for modeling
    required_cols = ['time', 'orgunitname', 'disease_cases', 'pop', 'temp_max', 'preci']
    missing = set(required_cols) - set(df.columns)
    if missing:
        print(f"Error: Missing required columns: {missing}")
        return
    
    print(f"Total rows in dataset: {len(df)}")
    
    # Drop current month (202511) as it has in-progress observations
    current_month = 202511
    df = df[df['time'] != current_month].copy()
    print(f"Rows after dropping current month ({current_month}): {len(df)}")
    
    # Create timeid (like in R script - based on sorted unique time values)
    sorted_times = sorted(df['time'].unique())
    time_to_id = {t: i+1 for i, t in enumerate(sorted_times)}
    df['timeid'] = df['time'].map(time_to_id)
    
    # Filter for timeid > 4 (to allow for lagged variables)
    df_filtered = df[df['timeid'] > 4].copy()
    
    print(f"Rows after filtering (timeid > 4): {len(df_filtered)}")
    
    # Compute lagged variables (needed for completeness check)
    print("Computing lagged climate variables...")
    
    def add_lags(data, col, max_lag):
        """Add lagged columns for a given variable."""
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
    
    # Center the climate variables (matching R code)
    df_filtered['PRCP'] = df_filtered['preci'] - df_filtered['preci'].mean()
    df_filtered['TEMPmax'] = df_filtered['temp_max'] - df_filtered['temp_max'].mean()
    
    # Add lags
    df_filtered = add_lags(df_filtered, 'PRCP', 2)
    df_filtered = add_lags(df_filtered, 'TEMPmax', 3)
    
    print("Lagged variables computed.")
    
    # Get all spatial units and time points
    all_units = sorted(df_filtered['orgunitname'].unique())
    all_times = sorted(df_filtered['timeid'].unique())
    
    id_to_time = {v: k for k, v in time_to_id.items()}
    most_recent_timeid = all_times[-1]
    most_recent_time = id_to_time[most_recent_timeid]
    
    # Define all variables needed for modeling (MUST be defined early)
    model_vars = ['disease_cases', 'pop', 'preci', 'temp_max',
                  'lag1_PRCP', 'lag2_PRCP', 
                  'lag1_TEMPmax', 'lag2_TEMPmax', 'lag3_TEMPmax']
    
    print(f"\nTotal spatial units: {len(all_units)}")
    print(f"Total time points (after filtering): {len(all_times)}")
    print(f"Time point range: {all_times[0]} to {all_times[-1]}")
    print(f"Most recent complete month: {most_recent_time} (timeid: {most_recent_timeid})")
    
    print(f"\nChecking completeness for {len(model_vars)} variables:")
    print(f"  {', '.join(model_vars)}")
    
    # For each time point, check if all units have complete data for ALL variables
    time_complete = []
    for t in all_times:
        time_data = df_filtered[df_filtered['timeid'] == t]
        
        # Check if all units are present
        units_present = set(time_data['orgunitname'].unique())
        all_present = units_present == set(all_units)
        
        # Check if ALL model variables are non-NA for all units
        if all_present:
            all_complete = time_data[model_vars].notna().all().all()
        else:
            all_complete = False
        
        time_complete.append(all_present and all_complete)
    
    time_complete = np.array(time_complete)
    
    print(f"\nTime points with complete data: {time_complete.sum()} out of {len(time_complete)}")
    
    # Analyze missing data patterns
    print("\nAnalyzing missing data patterns...")
    
    # Count how many units have complete data at each time point
    units_complete_per_time = []
    for t in all_times:
        time_data = df_filtered[df_filtered['timeid'] == t]
        units_present = set(time_data['orgunitname'].unique())
        cases_complete = time_data['disease_cases'].notna().sum()
        units_complete_per_time.append({
            'timeid': t,
            'units_present': len(units_present),
            'cases_complete': cases_complete
        })
    
    units_df = pd.DataFrame(units_complete_per_time)
    print(f"\nUnits with complete data per time point:")
    print(f"  Min: {units_df['cases_complete'].min()}")
    print(f"  Max: {units_df['cases_complete'].max()}")
    print(f"  Mean: {units_df['cases_complete'].mean():.1f}")
    
    # Find time point with most complete data
    best_time = units_df.loc[units_df['cases_complete'].idxmax()]
    print(f"\nBest time point: {best_time['timeid']} with {best_time['cases_complete']}/{len(all_units)} units having data")
    
    # Count data availability per unit
    print("\nAnalyzing per-unit data availability...")
    unit_completeness = []
    for unit in all_units:
        unit_data = df_filtered[df_filtered['orgunitname'] == unit]
        times_present = len(unit_data)
        # Check completeness for ALL model variables
        all_vars_complete = unit_data[model_vars].notna().all(axis=1).sum()
        unit_completeness.append({
            'orgunitname': unit,
            'times_present': times_present,
            'cases_complete': all_vars_complete
        })
    
    unit_df = pd.DataFrame(unit_completeness)
    print(f"\nData availability per unit:")
    print(f"  Units with all {len(all_times)} time points: {(unit_df['times_present'] == len(all_times)).sum()}")
    print(f"  Units with complete data (all variables) for all their time points: {(unit_df['cases_complete'] == unit_df['times_present']).sum()}")
    
    # Helper: find longest recent consecutive span (used for trade-off checks)
    def find_longest_recent_span_local(logical_array):
        """Find longest consecutive True span ending at the last index (or most recent True run)."""
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

    # --- Trade-off analysis: for many subset sizes, find longest recent span ---
    print("\nPerforming trade-off analysis (units vs. recent months)...")
    sorted_units = unit_df.sort_values('cases_complete', ascending=False)

    def check_subset_for_trade(n_units):
        """Return span dict for top n_units by coverage."""
        subset_units = set(sorted_units.head(n_units)['orgunitname'].tolist())
        subset_df = df_filtered[df_filtered['orgunitname'].isin(subset_units)]
        subset_complete = []
        for t in all_times:
            time_data = subset_df[subset_df['timeid'] == t]
            if len(time_data) != n_units:
                subset_complete.append(False)
                continue
            # Check ALL model variables are complete
            subset_complete.append(time_data[model_vars].notna().all().all())
        subset_complete = np.array(subset_complete)
        return find_longest_recent_span_local(subset_complete)

    trade_rows = []
    # We'll check every n from all_units down to 50 (or 1 if fewer units)
    min_units = max(1, int(len(all_units) * 0.3))
    for n in range(len(all_units), min_units - 1, -1):
        span = check_subset_for_trade(n)
        if span['length'] == 0:
            trade_rows.append({
                'n_units': n,
                'n_months': 0,
                'start_timeid': None,
                'end_timeid': None,
                'start_time': None,
                'end_time': None,
                'ends_at_recent': False,
                'total_obs': 0
            })
        else:
            stid = all_times[span['start']]
            enid = all_times[span['end']]
            st = id_to_time[stid]
            en = id_to_time[enid]
            ends_recent = (enid == most_recent_timeid)
            trade_rows.append({
                'n_units': n,
                'n_months': span['length'],
                'start_timeid': stid,
                'end_timeid': enid,
                'start_time': st,
                'end_time': en,
                'ends_at_recent': ends_recent,
                'total_obs': n * span['length']
            })

    trade_df = pd.DataFrame(trade_rows)
    trade_out = Path("output/tradeoff_orgunits_vs_months.csv")
    trade_out.parent.mkdir(exist_ok=True)
    trade_df.to_csv(trade_out, index=False)
    print(f"Saved trade-off results to: {trade_out}")

    # Find units with least data
    worst_units = unit_df.nsmallest(5, 'cases_complete')
    print(f"\nUnits with least complete disease_cases data:")
    for _, row in worst_units.iterrows():
        print(f"  {row['orgunitname']}: {row['cases_complete']}/{row['times_present']} time points")
    
    # Find longest consecutive span ending at most recent time
    # Search backwards from the most recent time point
    print("\nSearching for longest recent consecutive span (moving backwards from most recent month)...")
    
    def find_longest_recent_span(logical_array):
        """Find longest consecutive True span ending at the last index."""
        if not logical_array.any():
            return {'start': None, 'end': None, 'length': 0}
        
        # Search backwards from the end
        end_idx = len(logical_array) - 1
        
        # Find the longest consecutive True sequence ending at end_idx
        length = 0
        for i in range(end_idx, -1, -1):
            if logical_array[i]:
                length += 1
            else:
                break
        
        if length > 0:
            return {'start': end_idx - length + 1, 'end': end_idx, 'length': length}
        
        # If the most recent time doesn't have complete data, find the most recent complete span
        # Find the last True value
        true_indices = np.where(logical_array)[0]
        if len(true_indices) == 0:
            return {'start': None, 'end': None, 'length': 0}
        
        # Work backwards from the last True value
        last_true = true_indices[-1]
        length = 1
        for i in range(last_true - 1, -1, -1):
            if logical_array[i]:
                length += 1
            else:
                break
        
        return {'start': last_true - length + 1, 'end': last_true, 'length': length}
    
    span_info = find_longest_recent_span(time_complete)
    
    if span_info['length'] == 0:
        print("\n" + "="*60)
        print("RESULT: No time span found where ALL spatial units have complete disease_cases.")
        print("="*60)
        
        # Find the largest subset of units that have a complete time span
        print("\nSearching for largest subset of units with complete time span...")
        
        # Start by excluding units with worst data coverage
        sorted_units = unit_df.sort_values('cases_complete', ascending=False)
        
        # Binary search for the largest subset
        left, right = 1, len(all_units)
        best_result = None
        
        def check_subset(n_units):
            """Check if n_units with best coverage have a complete recent time span."""
            subset_units = set(sorted_units.head(n_units)['orgunitname'].tolist())
            
            # Create a matrix: time x unit x has_data
            subset_df = df_filtered[df_filtered['orgunitname'].isin(subset_units)]
            
            subset_complete = []
            for t in all_times:
                time_data = subset_df[subset_df['timeid'] == t]
                
                # Quick check: do we have n_units entries?
                if len(time_data) != n_units:
                    subset_complete.append(False)
                    continue
                
                # Check all have complete data for ALL model variables
                all_complete = time_data[model_vars].notna().all().all()
                subset_complete.append(all_complete)
            
            subset_complete = np.array(subset_complete)
            return find_longest_recent_span(subset_complete)
        
        print("  Searching across different subset sizes...")
        
        # Try different subset sizes and track all results
        all_results = []
        
        # Check several key subset sizes
        test_sizes = []
        for pct in [100, 99, 98, 95, 90, 85, 80, 75, 70]:
            size = int(len(all_units) * pct / 100)
            if size > 0 and size not in test_sizes:
                test_sizes.append(size)
        
        for n_units in test_sizes:
            span = check_subset(n_units)
            if span['length'] > 0:
                all_results.append({
                    'n_units': n_units,
                    'units': sorted_units.head(n_units)['orgunitname'].tolist(),
                    'span': span,
                    'times': all_times
                })
        
        if not all_results:
            # Binary search for any valid result
            while left <= right:
                mid = (left + right) // 2
                span = check_subset(mid)
                
                if span['length'] > 0:
                    all_results.append({
                        'n_units': mid,
                        'units': sorted_units.head(mid)['orgunitname'].tolist(),
                        'span': span,
                        'times': all_times
                    })
                    left = mid + 1
                else:
                    right = mid - 1
        
        # Find best result prioritizing: 1) ends at recent, 2) longer span, 3) more units
        if all_results:
            def score_result(r):
                end_timeid = r['times'][r['span']['end']]
                is_recent = (end_timeid == most_recent_timeid)
                return (
                    1000000 if is_recent else 0,  # Strongly prefer recent
                    r['span']['length'],            # Then prefer longer
                    r['n_units']                    # Then prefer more units
                )
            
            all_results.sort(key=score_result, reverse=True)
            best_result = all_results[0]
            
            # Also show top 3 alternatives
            print(f"\n  Found {len(all_results)} valid configurations")
            print(f"  Top 3 options:")
            for i, r in enumerate(all_results[:3], 1):
                end_id = r['times'][r['span']['end']]
                start_id = r['times'][r['span']['start']]
                is_recent = (end_id == most_recent_timeid)
                months_behind = most_recent_timeid - end_id
                print(f"    {i}. {r['n_units']:3d} units × {r['span']['length']:3d} months = {r['n_units'] * r['span']['length']:5d} obs "
                      f"(timeid {start_id}-{end_id}, {'RECENT' if is_recent else f'{months_behind} mo. old'})")
        else:
            best_result = None
        
        if best_result:
            start_timeid = best_result['times'][best_result['span']['start']]
            end_timeid = best_result['times'][best_result['span']['end']]
            start_time = id_to_time[start_timeid]
            end_time = id_to_time[end_timeid]
            
            is_ending_at_recent = (end_timeid == most_recent_timeid)
            
            print("\n" + "="*60)
            print("ALTERNATIVE: LARGEST SUBSET WITH COMPLETE RECENT TIME SPAN")
            print("="*60)
            print(f"Number of units: {best_result['n_units']} (out of {len(all_units)})")
            print(f"Start time ID:   {start_timeid} (time: {start_time})")
            print(f"End time ID:     {end_timeid} (time: {end_time})")
            print(f"Length:          {best_result['span']['length']} time points")
            print(f"Total obs:       {best_result['n_units']} × {best_result['span']['length']} = {best_result['n_units'] * best_result['span']['length']}")
            if is_ending_at_recent:
                print(f"Status:          ✓ Ends at most recent month ({end_time})")
            else:
                print(f"Status:          ⚠ Ends {most_recent_timeid - end_timeid} months before most recent")
            print("="*60)
            
            # Show which units are excluded
            excluded_units = sorted(set(all_units) - set(best_result['units']))
            print(f"\nExcluded {len(excluded_units)} units:")
            for unit in excluded_units:
                unit_info = unit_df[unit_df['orgunitname'] == unit].iloc[0]
                print(f"  - {unit}: {unit_info['cases_complete']}/{unit_info['times_present']} complete")
            
            # Save recommended units to CSV
            output_path = Path("output/recommended_orgunits.csv")
            output_path.parent.mkdir(exist_ok=True)
            
            recommended_df = pd.DataFrame({
                'orgunitname': best_result['units'],
                'start_time': start_time,
                'end_time': end_time,
                'start_timeid': start_timeid,
                'end_timeid': end_timeid,
                'n_months': best_result['span']['length']
            })
            recommended_df.to_csv(output_path, index=False)
            print(f"\n✓ Saved {len(best_result['units'])} recommended orgunits to: {output_path}")
        else:
            print("\nCould not find any subset of units with a complete time span.")
        
        print("\nOther suggestions:")
        print("1. Impute missing values for disease_cases")
        print("2. Use a different approach that handles missing data")
        return
    
    # Get actual time IDs for this span
    start_timeid = all_times[span_info['start']]
    end_timeid = all_times[span_info['end']]
    
    # Get corresponding time values (e.g., 201701, 201702, etc.)
    start_time = id_to_time[start_timeid]
    end_time = id_to_time[end_timeid]
    
    is_ending_at_recent = (end_timeid == most_recent_timeid)
    
    print("\n" + "="*60)
    print("LONGEST COMPLETE RECENT TIME SPAN FOUND (ALL UNITS):")
    print("="*60)
    print(f"Start time ID: {start_timeid} (time: {start_time})")
    print(f"End time ID:   {end_timeid} (time: {end_time})")
    print(f"Length:        {span_info['length']} time points")
    if is_ending_at_recent:
        print(f"Status:        ✓ Ends at most recent month ({end_time})")
    else:
        print(f"Status:        ⚠ Ends {most_recent_timeid - end_timeid} months before most recent")
    print("="*60)
    
    # Count total observations in this span
    span_data = df_filtered[
        (df_filtered['timeid'] >= start_timeid) & 
        (df_filtered['timeid'] <= end_timeid)
    ]
    print(f"\nTotal observations in this span: {len(span_data)}")
    print(f"Expected observations: {len(all_units)} units × {span_info['length']} times = {len(all_units) * span_info['length']}")
    
    # Show which time points are complete (for verification)
    print("\nTime points with complete data for all units:")
    complete_times = [all_times[i] for i, complete in enumerate(time_complete) if complete]
    if len(complete_times) <= 20:
        print(f"  {complete_times}")
    else:
        print(f"  First 10: {complete_times[:10]}")
        print(f"  Last 10:  {complete_times[-10:]}")
    
    # Save recommended units to CSV
    output_path = Path("output/recommended_orgunits.csv")
    output_path.parent.mkdir(exist_ok=True)
    
    recommended_df = pd.DataFrame({
        'orgunitname': all_units,
        'start_time': start_time,
        'end_time': end_time,
        'start_timeid': start_timeid,
        'end_timeid': end_timeid,
        'n_months': span_info['length']
    })
    recommended_df.to_csv(output_path, index=False)
    print(f"\n✓ Saved {len(all_units)} recommended orgunits to: {output_path}")


if __name__ == "__main__":
    main()
