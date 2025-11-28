#!/usr/bin/env python3
"""
DBH Comparison Analysis
Compares different DBH measurement methods:
- 3DFIN DBH
- Field D (reference)
- TLS DBH
- DBH (Point Cloud Tools)
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats

# Set style for better-looking plots
sns.set_style("whitegrid")
plt.rcParams['figure.figsize'] = (16, 12)

# Load data
print("Loading data...")
df = pd.read_csv('matlab_treelist_combined_with_threedfin_A.csv')

# Convert threedfin_dbh from meters to cm (if needed)
# Based on the MATLAB script, threedfin_dbh is multiplied by 100
df['threedfin_dbh_cm'] = df['threedfin_dbh'] * 100

# Convert TLS_dbh from meters to cm
df['TLS_dbh_cm'] = df['TLS_dbh'] * 100

# Create a clean dataset with all four methods
methods_df = df[['field_D', 'TLS_dbh_cm', 'dbh', 'threedfin_dbh_cm']].copy()
methods_df.columns = ['Field_D', 'TLS_DBH', 'PC_Tools_DBH', '3DFIN_DBH']

# Remove rows with all NaN
methods_df = methods_df.dropna(how='all')

print(f"\nTotal records: {len(methods_df)}")
print(f"\nData availability:")
print(methods_df.count())

# Create figure with subplots
fig = plt.figure(figsize=(20, 16))

# 1. Scatter plots comparing each method to Field_D (reference)
methods_to_compare = ['TLS_DBH', 'PC_Tools_DBH', '3DFIN_DBH']
colors = ['#2E86AB', '#A23B72', '#F18F01']

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
        ax.set_title(f'{method.replace("_", " ")} vs Field D', fontsize=14, fontweight='bold')
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
        ax.set_title(f'{method.replace("_", " ")} Residuals', fontsize=14, fontweight='bold')
        ax.legend()
        ax.grid(True, alpha=0.3)

# 3. Distribution comparison
ax = plt.subplot(3, 3, 7)
for method, color in zip(['Field_D'] + methods_to_compare, ['black'] + colors):
    valid_values = methods_df[method].dropna()
    if len(valid_values) > 0:
        ax.hist(valid_values, bins=30, alpha=0.4, label=method.replace('_', ' '),
                color=color, edgecolor='black', linewidth=1)

ax.set_xlabel('DBH (cm)', fontsize=12, fontweight='bold')
ax.set_ylabel('Frequency', fontsize=12, fontweight='bold')
ax.set_title('DBH Distribution Comparison', fontsize=14, fontweight='bold')
ax.legend()
ax.grid(True, alpha=0.3)

# 4. Box plot comparison
ax = plt.subplot(3, 3, 8)
box_data = [methods_df[col].dropna() for col in methods_df.columns]
bp = ax.boxplot(box_data, tick_labels=[col.replace('_', '\n') for col in methods_df.columns],
                patch_artist=True)

# Color the boxes
for patch, color in zip(bp['boxes'], ['black'] + colors):
    patch.set_facecolor(color)
    patch.set_alpha(0.6)

ax.set_ylabel('DBH (cm)', fontsize=12, fontweight='bold')
ax.set_title('DBH Distribution by Method', fontsize=14, fontweight='bold')
ax.grid(True, alpha=0.3, axis='y')

# 5. Correlation matrix
ax = plt.subplot(3, 3, 9)
corr_matrix = methods_df.corr()
sns.heatmap(corr_matrix, annot=True, fmt='.3f', cmap='coolwarm',
            square=True, cbar_kws={'label': 'Correlation'}, ax=ax,
            vmin=0, vmax=1, linewidths=1, linecolor='black')
ax.set_title('Correlation Matrix', fontsize=14, fontweight='bold')

plt.tight_layout()
plt.savefig('dbh_comparison_analysis.png', dpi=300, bbox_inches='tight')
print("\n✓ Saved: dbh_comparison_analysis.png")

# Create a second figure for pairwise comparisons
fig2, axes = plt.subplots(4, 4, figsize=(20, 20))

all_methods = ['Field_D', 'TLS_DBH', 'PC_Tools_DBH', '3DFIN_DBH']
for i, method1 in enumerate(all_methods):
    for j, method2 in enumerate(all_methods):
        ax = axes[i, j]

        if i == j:
            # Diagonal: show distribution
            valid_values = methods_df[method1].dropna()
            if len(valid_values) > 0:
                ax.hist(valid_values, bins=30, alpha=0.7, color=colors[min(i, len(colors)-1)],
                       edgecolor='black', linewidth=1)
                ax.set_title(f'{method1.replace("_", " ")}', fontweight='bold')
                ax.set_ylabel('Frequency')
        else:
            # Off-diagonal: scatter plot
            valid_data = methods_df[[method1, method2]].dropna()
            if len(valid_data) > 0:
                x = valid_data[method1]
                y = valid_data[method2]

                ax.scatter(x, y, alpha=0.4, s=20, edgecolors='black', linewidth=0.5)

                # 1:1 line
                max_val = max(x.max(), y.max())
                min_val = min(x.min(), y.min())
                ax.plot([min_val, max_val], [min_val, max_val], 'r--', linewidth=1.5, alpha=0.7)

                # Correlation
                r_value = methods_df[[method1, method2]].corr().iloc[0, 1]
                ax.text(0.05, 0.95, f'r = {r_value:.3f}', transform=ax.transAxes,
                       verticalalignment='top', bbox=dict(boxstyle='round', facecolor='yellow', alpha=0.7))

                if i == 3:
                    ax.set_xlabel(method1.replace('_', ' '), fontweight='bold')
                if j == 0:
                    ax.set_ylabel(method2.replace('_', ' '), fontweight='bold')

        ax.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('dbh_pairwise_comparison.png', dpi=300, bbox_inches='tight')
print("✓ Saved: dbh_pairwise_comparison.png")

# Create summary statistics table
print("\n" + "="*80)
print("SUMMARY STATISTICS")
print("="*80)
print("\nDescriptive Statistics (cm):")
print(methods_df.describe().round(2))

print("\n" + "-"*80)
print("Correlation with Field_D (reference):")
print("-"*80)
for method in methods_to_compare:
    valid_data = methods_df[['Field_D', method]].dropna()
    if len(valid_data) > 0:
        x = valid_data['Field_D']
        y = valid_data[method]

        r_value = np.corrcoef(x, y)[0, 1]
        rmse = np.sqrt(np.mean((x - y)**2))
        bias = np.mean(y - x)
        mae = np.mean(np.abs(x - y))

        print(f"\n{method.replace('_', ' ')}:")
        print(f"  Sample size: {len(valid_data)}")
        print(f"  Correlation (r): {r_value:.4f}")
        print(f"  R²: {r_value**2:.4f}")
        print(f"  RMSE: {rmse:.3f} cm")
        print(f"  Bias: {bias:.3f} cm")
        print(f"  MAE: {mae:.3f} cm")

# Save statistics to CSV
stats_summary = []
for method in methods_to_compare:
    valid_data = methods_df[['Field_D', method]].dropna()
    if len(valid_data) > 0:
        x = valid_data['Field_D']
        y = valid_data[method]

        r_value = np.corrcoef(x, y)[0, 1]
        rmse = np.sqrt(np.mean((x - y)**2))
        bias = np.mean(y - x)
        mae = np.mean(np.abs(x - y))

        stats_summary.append({
            'Method': method,
            'N': len(valid_data),
            'Correlation_r': r_value,
            'R_squared': r_value**2,
            'RMSE_cm': rmse,
            'Bias_cm': bias,
            'MAE_cm': mae,
            'Mean_cm': y.mean(),
            'Std_cm': y.std()
        })

stats_df = pd.DataFrame(stats_summary)
stats_df.to_csv('dbh_comparison_statistics.csv', index=False)
print("\n✓ Saved: dbh_comparison_statistics.csv")

print("\n" + "="*80)
print("Analysis complete! Generated files:")
print("  - dbh_comparison_analysis.png")
print("  - dbh_pairwise_comparison.png")
print("  - dbh_comparison_statistics.csv")
print("="*80)
