#!/usr/bin/env python3
"""
DBH Comparison Analysis - Enhanced MATLAB-style version
Compares different DBH measurement methods following MATLAB workflow:
- Multiple trajectories (A=loop, L=line)
- 3DFIN DBH vs Field D (reference)
- TLS DBH vs Field D
- DBH (Point Cloud Tools) vs Field D
- Multiple threshold analysis (>=5, >=7, >=10 cm)
- Plot-level aggregation
- Tree finding ratios (detection rates)
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
import os

# Set style for better-looking plots
sns.set_style("whitegrid")
plt.rcParams['figure.figsize'] = (20, 16)

#%% Configuration
# Trajectories to analyze (like MATLAB: traj = {'A', 'L'})
trajectories = []
if os.path.exists('matlab_treelist_combined_with_threedfin_A.csv'):
    trajectories.append('A')  # Loop trajectory
if os.path.exists('matlab_treelist_combined_with_threedfin_L.csv'):
    trajectories.append('L')  # Line trajectory

if not trajectories:
    print("ERROR: No trajectory files found!")
    exit(1)

print(f"Found trajectories: {trajectories}")

# Analysis parameters (following MATLAB style)
distance_threshold = 7.5  # cm - IMPORTANT: Distance correlates with accuracy
thresholds = [5, 7, 10]  # DBH thresholds for detailed analysis (cm)

# Storage structures (like MATLAB structures)
whole_population = {}
plot_levels = {}
stats_plot_levels = {}
correlation_tables = {}
All_measured_trees = {}

#%% Loop over trajectories (following MATLAB: for i = 1:length(traj))
for walk in trajectories:
    print("\n" + "="*80)
    print(f"Processing trajectory: {walk}")
    print("="*80)

    #%% Load data
    filename = f'matlab_treelist_combined_with_threedfin_{walk}.csv'
    print(f"\nLoading: {filename}")
    df = pd.read_csv(filename)

    # Convert threedfin_dbh from meters to cm
    # Following MATLAB: matlab_treelist_combined_all_total.threedfin_dbh = matlab_treelist_combined_all_total.threedfin_dbh * 100;
    df['threedfin_dbh_cm'] = df['threedfin_dbh'] * 100

    # Convert TLS_dbh from meters to cm
    df['TLS_dbh_cm'] = df['TLS_dbh'] * 100

    #%% IMPORTANT: Filter 3DFIN data based on distance threshold
    # Following MATLAB: matlab_treelist_combined_all_total = matlab_treelist_combined_all_total(matlab_treelist_combined_all_total.distance_to_threedfin_cm < 7.5,:);
    # Distance correlates with measurement accuracy - only keep close matches
    print(f"\nOriginal 3DFIN records: {df['threedfin_dbh_cm'].notna().sum()}")
    threedfin_filtered = df['distance_to_threedfin_cm'] < distance_threshold

    # Store unfiltered version for detection rate calculation
    df_unfiltered = df.copy()

    # Apply distance filter to 3DFIN data only
    df.loc[~threedfin_filtered, 'threedfin_dbh_cm'] = np.nan
    print(f"3DFIN records after distance filter (<{distance_threshold} cm): {df['threedfin_dbh_cm'].notna().sum()}")
    print(f"Filtered out {(~threedfin_filtered & df['distance_to_threedfin_cm'].notna()).sum()} trees with distance >= {distance_threshold} cm")

    # Filter based on distance (following MATLAB)
    valid_idx = ~df['threedfin_dbh_cm'].isna() & (df['threedfin_dbh_cm'] != 0)
    df_measured = df[valid_idx].copy()

    #%% Calculate tree finding ratios // Dr // Detection rate against reference
    # Following MATLAB calculations
    total_trees_amount = len(df)  # For detection rate (would be from field reference in full version)

    # Trees which have matchID
    threedfin_found_trees = df_unfiltered['threedfin_match'].notna().sum()
    tree_segmented_ratio = (threedfin_found_trees / total_trees_amount) * 100

    # Trees with valid DBH measurements
    threedfin_measured_trees = ((~df_unfiltered['threedfin_dbh_cm'].isna()) &
                                (df_unfiltered['threedfin_dbh_cm'] != 0)).sum()
    threedfin_measured_ratio = (threedfin_measured_trees / total_trees_amount) * 100

    # Trees after filtering
    threedfin_total_trees_amount = len(df_measured)
    tree_finding_ratio = (threedfin_total_trees_amount / total_trees_amount) * 100

    print(f"\n{walk} Tree Finding Ratio: {tree_finding_ratio:.2f}%")
    print(f"{walk} Tree Segmentation Ratio: {tree_segmented_ratio:.2f}%")
    print(f"{walk} Tree Measured Ratio: {threedfin_measured_ratio:.2f}%")

    #%% Create clean dataset with all four methods
    methods_df = df[['field_D', 'TLS_dbh_cm', 'dbh', 'threedfin_dbh_cm', 'field_KoealaId']].copy()
    methods_df.columns = ['Field_D', 'TLS_DBH', 'PC_Tools_DBH', '3DFIN_DBH', 'PlotID']

    # Remove rows with all NaN
    methods_df = methods_df.dropna(how='all')

    print(f"\nTotal records: {len(methods_df)}")
    print(f"\nData availability:")
    print(methods_df.drop('PlotID', axis=1).count())

    #%% DBH: RMSE and Bias for DBH at thresholds >= 5, 7, and 10 cm for whole population
    # Following MATLAB structure exactly

    methods_to_compare = ['TLS_DBH', 'PC_Tools_DBH', '3DFIN_DBH']

    # Initialize results structure (following MATLAB: DBH_results = struct(...))
    DBH_results = {method: [] for method in methods_to_compare}

    print(f"\n{'='*80}")
    print(f"THRESHOLD ANALYSIS - Trajectory {walk}")
    print(f"{'='*80}")

    for method in methods_to_compare:
        print(f"\n{method.replace('_', ' ')}:")
        print(f"{'Threshold':<12} {'N':<8} {'RMSE (cm)':<12} {'Bias (cm)':<12} {'R²':<10} {'r':<10}")
        print("-" * 70)

        for t in thresholds:
            # Following MATLAB mask creation
            mask = (methods_df['Field_D'] >= t) & \
                   (methods_df[method] > 0) & \
                   (methods_df['Field_D'] > 0) & \
                   (~methods_df['Field_D'].isna()) & \
                   (~methods_df[method].isna())

            fieldDBH = methods_df.loc[mask, 'Field_D']
            methodDBH = methods_df.loc[mask, method]

            if len(fieldDBH) > 0:
                # Calculate RMSE and Bias (following MATLAB exactly)
                RMSE_dbh = np.sqrt(np.mean((fieldDBH - methodDBH)**2))
                Bias_dbh = np.mean(methodDBH - fieldDBH)
                mean_field_D = np.mean(fieldDBH)
                RMSE_dbh_percent = (RMSE_dbh / mean_field_D) * 100

                # Calculate correlation
                r_dbh = np.corrcoef(methodDBH, fieldDBH)[0, 1]
                r_dbh_2 = r_dbh**2

                # Calculate R² (coefficient of determination)
                SS_res = np.sum((fieldDBH - methodDBH)**2)
                SS_tot = np.sum((fieldDBH - np.mean(fieldDBH))**2)
                r2_dbh = 1 - (SS_res / SS_tot)

                # Store results
                DBH_results[method].append({
                    'Threshold': t,
                    'RMSE': RMSE_dbh,
                    'Bias': Bias_dbh,
                    'SampleSize': len(fieldDBH),
                    'RMSE_percent': RMSE_dbh_percent,
                    'r': r_dbh,
                    'r2': r_dbh_2,
                    'R2': r2_dbh
                })

                print(f">= {t} cm    {len(fieldDBH):<8} {RMSE_dbh:<12.3f} {Bias_dbh:<12.3f} {r2_dbh:<10.4f} {r_dbh:<10.4f}")

    # Store in whole_population structure (following MATLAB)
    whole_population[walk] = {'DBH': DBH_results}

    #%% Plot-level aggregation (following MATLAB: groupsummary)
    print(f"\n{'='*80}")
    print(f"PLOT-LEVEL ANALYSIS - Trajectory {walk}")
    print(f"{'='*80}")

    # Get unique plots
    unique_plots = methods_df['PlotID'].dropna().unique()
    print(f"\nNumber of plots: {len(unique_plots)}")

    # Initialize plot-level results
    plot_level_results = []

    for plot_id in unique_plots:
        plot_data = methods_df[methods_df['PlotID'] == plot_id]

        result = {'PlotID': plot_id}

        # Calculate statistics for each method (following MATLAB groupsummary)
        for method in ['Field_D', 'TLS_DBH', 'PC_Tools_DBH', '3DFIN_DBH']:
            valid_values = plot_data[method].dropna()

            if len(valid_values) > 0:
                result[f'{method}_mean'] = valid_values.mean()
                result[f'{method}_std'] = valid_values.std()
                result[f'{method}_min'] = valid_values.min()
                result[f'{method}_max'] = valid_values.max()
                result[f'{method}_count'] = len(valid_values)
            else:
                result[f'{method}_mean'] = np.nan
                result[f'{method}_std'] = np.nan
                result[f'{method}_min'] = np.nan
                result[f'{method}_max'] = np.nan
                result[f'{method}_count'] = 0

        # Calculate RMSE and Bias per plot for each method
        for method in methods_to_compare:
            valid_idx = ~plot_data['Field_D'].isna() & ~plot_data[method].isna()
            if valid_idx.sum() > 0:
                field_vals = plot_data.loc[valid_idx, 'Field_D']
                method_vals = plot_data.loc[valid_idx, method]

                rmse = np.sqrt(np.mean((field_vals - method_vals)**2))
                bias = np.mean(method_vals - field_vals)

                result[f'{method}_RMSE'] = rmse
                result[f'{method}_Bias'] = bias
            else:
                result[f'{method}_RMSE'] = np.nan
                result[f'{method}_Bias'] = np.nan

        plot_level_results.append(result)

    # Convert to DataFrame (like MATLAB table)
    plot_level = pd.DataFrame(plot_level_results)
    plot_levels[walk] = plot_level

    # Calculate overall plot-level statistics
    stats_plot_level = {}
    for method in methods_to_compare:
        valid_plots = plot_level[[f'Field_D_mean', f'{method}_mean']].dropna()

        if len(valid_plots) > 0:
            field_means = valid_plots['Field_D_mean']
            method_means = valid_plots[f'{method}_mean']

            stats_plot_level[f'{method}_RMSE'] = np.sqrt(np.mean((field_means - method_means)**2))
            stats_plot_level[f'{method}_Bias'] = np.mean(method_means - field_means)
            stats_plot_level[f'{method}_r'] = np.corrcoef(method_means, field_means)[0, 1]
            stats_plot_level[f'{method}_r2'] = stats_plot_level[f'{method}_r']**2

    stats_plot_levels[walk] = stats_plot_level

    print(f"\nPlot-level correlations:")
    for method in methods_to_compare:
        if f'{method}_r' in stats_plot_level:
            print(f"  {method.replace('_', ' ')}: r = {stats_plot_level[f'{method}_r']:.4f}, R² = {stats_plot_level[f'{method}_r2']:.4f}")

    #%% Store all measured trees (following MATLAB)
    All_measured_trees[walk] = df_measured

    #%% Create visualization for this trajectory
    fig = plt.figure(figsize=(20, 16))
    colors = ['#2E86AB', '#A23B72', '#F18F01']

    # 1. Scatter plots comparing each method to Field_D (reference)
    for idx, (method, color) in enumerate(zip(methods_to_compare, colors), 1):
        ax = plt.subplot(3, 3, idx)

        # Filter valid data
        valid_data = methods_df[['Field_D', method]].dropna()

        if len(valid_data) > 0:
            x = valid_data['Field_D']
            y = valid_data[method]

            # Scatter plot
            ax.scatter(x, y, alpha=0.5, s=30, color=color, edgecolors='black', linewidth=0.5)

            # 1:1 line
            max_val = max(x.max(), y.max())
            min_val = min(x.min(), y.min())
            ax.plot([min_val, max_val], [min_val, max_val], 'k--', linewidth=2, label='1:1 line')

            # Linear regression
            slope, intercept, r_value, p_value, std_err = stats.linregress(x, y)
            line_x = np.array([min_val, max_val])
            line_y = slope * line_x + intercept
            ax.plot(line_x, line_y, 'r-', linewidth=2, label=f'Fit: y={slope:.2f}x+{intercept:.2f}')

            # Calculate RMSE and Bias
            rmse = np.sqrt(np.mean((x - y)**2))
            bias = np.mean(y - x)

            # Add statistics
            stats_text = f'n = {len(valid_data)}\n'
            stats_text += f'R² = {r_value**2:.3f}\n'
            stats_text += f'RMSE = {rmse:.2f} cm\n'
            stats_text += f'Bias = {bias:.2f} cm'

            ax.text(0.05, 0.95, stats_text, transform=ax.transAxes,
                    verticalalignment='top', bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8),
                    fontsize=10)

            ax.set_xlabel('Field D (cm)', fontsize=12, fontweight='bold')
            ax.set_ylabel(f'{method.replace("_", " ")} (cm)', fontsize=12, fontweight='bold')

            # Add distance filter note for 3DFIN
            title = f'{method.replace("_", " ")} vs Field D - Walk {walk}'
            if method == '3DFIN_DBH':
                title += f'\n(distance < {distance_threshold} cm)'
            ax.set_title(title, fontsize=14, fontweight='bold')

            ax.legend(loc='lower right')
            ax.grid(True, alpha=0.3)

    # 2. Residual plots
    for idx, (method, color) in enumerate(zip(methods_to_compare, colors), 4):
        ax = plt.subplot(3, 3, idx)

        valid_data = methods_df[['Field_D', method]].dropna()

        if len(valid_data) > 0:
            x = valid_data['Field_D']
            residuals = valid_data[method] - x

            ax.scatter(x, residuals, alpha=0.5, s=30, color=color, edgecolors='black', linewidth=0.5)
            ax.axhline(y=0, color='k', linestyle='--', linewidth=2)

            # Add mean residual line
            mean_residual = residuals.mean()
            ax.axhline(y=mean_residual, color='r', linestyle='-', linewidth=2,
                       label=f'Mean bias = {mean_residual:.2f} cm')

            ax.set_xlabel('Field D (cm)', fontsize=12, fontweight='bold')
            ax.set_ylabel('Residual (cm)', fontsize=12, fontweight='bold')

            # Add distance filter note for 3DFIN
            title = f'{method.replace("_", " ")} Residuals - Walk {walk}'
            if method == '3DFIN_DBH':
                title += f'\n(distance < {distance_threshold} cm)'
            ax.set_title(title, fontsize=14, fontweight='bold')

            ax.legend()
            ax.grid(True, alpha=0.3)

    # 3. Distribution comparison
    ax = plt.subplot(3, 3, 7)
    for method, color in zip(['Field_D'] + methods_to_compare, ['black'] + colors):
        valid_values = methods_df[method].dropna()
        if len(valid_values) > 0:
            label = method.replace('_', ' ')
            if method == '3DFIN_DBH':
                label += f' (dist<{distance_threshold}cm)'
            ax.hist(valid_values, bins=30, alpha=0.4, label=label,
                    color=color, edgecolor='black', linewidth=1)

    ax.set_xlabel('DBH (cm)', fontsize=12, fontweight='bold')
    ax.set_ylabel('Frequency', fontsize=12, fontweight='bold')
    ax.set_title(f'DBH Distribution Comparison - Walk {walk}', fontsize=14, fontweight='bold')
    ax.legend()
    ax.grid(True, alpha=0.3)

    # 4. Box plot comparison
    ax = plt.subplot(3, 3, 8)
    box_data = [methods_df[col].dropna() for col in ['Field_D', 'TLS_DBH', 'PC_Tools_DBH', '3DFIN_DBH']]
    box_labels = []
    for col in ['Field_D', 'TLS_DBH', 'PC_Tools_DBH', '3DFIN_DBH']:
        label = col.replace('_', '\n')
        if col == '3DFIN_DBH':
            label += f'\n(d<{distance_threshold})'
        box_labels.append(label)
    bp = ax.boxplot(box_data, tick_labels=box_labels, patch_artist=True)

    # Color the boxes
    for patch, color in zip(bp['boxes'], ['black'] + colors):
        patch.set_facecolor(color)
        patch.set_alpha(0.6)

    ax.set_ylabel('DBH (cm)', fontsize=12, fontweight='bold')
    ax.set_title(f'DBH Distribution by Method - Walk {walk}', fontsize=14, fontweight='bold')
    ax.grid(True, alpha=0.3, axis='y')

    # 5. Correlation matrix
    ax = plt.subplot(3, 3, 9)
    corr_data = methods_df[['Field_D', 'TLS_DBH', 'PC_Tools_DBH', '3DFIN_DBH']].copy()
    corr_matrix = corr_data.corr()
    sns.heatmap(corr_matrix, annot=True, fmt='.3f', cmap='coolwarm',
                square=True, cbar_kws={'label': 'Correlation'}, ax=ax,
                vmin=0, vmax=1, linewidths=1, linecolor='black')
    ax.set_title(f'Correlation Matrix - Walk {walk}', fontsize=14, fontweight='bold')

    plt.tight_layout()
    plt.savefig(f'dbh_comparison_analysis_{walk}.png', dpi=300, bbox_inches='tight')
    print(f"\n✓ Saved: dbh_comparison_analysis_{walk}.png")

#%% Summary across all trajectories
print("\n" + "="*80)
print("SUMMARY ACROSS ALL TRAJECTORIES")
print("="*80)

# Create summary comparison table
summary_data = []
for walk in trajectories:
    if walk in whole_population:
        for method in methods_to_compare:
            for threshold_data in whole_population[walk]['DBH'][method]:
                summary_data.append({
                    'Trajectory': walk,
                    'Method': method,
                    'Threshold': threshold_data['Threshold'],
                    'N': threshold_data['SampleSize'],
                    'RMSE_cm': threshold_data['RMSE'],
                    'Bias_cm': threshold_data['Bias'],
                    'RMSE_percent': threshold_data['RMSE_percent'],
                    'r': threshold_data['r'],
                    'R2': threshold_data['R2']
                })

summary_df = pd.DataFrame(summary_data)
summary_df.to_csv('dbh_comparison_summary_all_trajectories.csv', index=False)
print("\n✓ Saved: dbh_comparison_summary_all_trajectories.csv")

# Export plot-level results
for walk in trajectories:
    if walk in plot_levels:
        plot_levels[walk].to_csv(f'plot_level_results_{walk}.csv', index=False)
        print(f"✓ Saved: plot_level_results_{walk}.csv")

print("\n" + "="*80)
print("Analysis complete!")
print("="*80)
print("\nGenerated files:")
for walk in trajectories:
    print(f"  - dbh_comparison_analysis_{walk}.png")
    print(f"  - plot_level_results_{walk}.csv")
print(f"  - dbh_comparison_summary_all_trajectories.csv")
