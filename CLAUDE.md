# CLAUDE.md - AI Assistant Guide for silver-bassoon

## Project Overview

**Project Name:** silver-bassoon
**Repository:** https://github.com/LLiiko/silver-bassoon.git
**Type:** Forestry Research & Data Analysis
**Primary Language:** MATLAB
**Domain:** Terrestrial Laser Scanning (TLS) & Tree Detection Analysis

### Purpose

This repository contains research code for analyzing tree measurements from Terrestrial Laser Scanning (TLS) data combined with 3DFIN (tree detection software) results. The analysis compares automated measurements against field reference data for forestry plot assessments.

## Repository Structure

```
silver-bassoon/
├── Train                                              # Repository URL reference
├── matlab_treelist_combined_with_threedfin_A.csv     # Primary dataset (~1.75MB)
└── threedfin_results_A.m                             # Main analysis script (~40KB)
```

### File Descriptions

#### 1. `Train`
- Simple text file containing the GitHub repository URL
- No functional code, serves as metadata

#### 2. `matlab_treelist_combined_with_threedfin_A.csv`
- **Size:** ~1.75 MB (1,749,721 bytes)
- **Format:** CSV with 47 columns
- **Content:** Tree-level measurements combining:
  - TLS (Terrestrial Laser Scanning) measurements
  - Field reference measurements
  - 3DFIN detection results

**Key Columns:**
- `TLS_tree_id`, `TLS_x`, `TLS_y`, `TLS_z` - Tree identification and coordinates
- `TLS_dbh`, `TLS_tree_ht` - TLS-measured diameter and height
- `field_D`, `field_H` - Field-measured diameter and height
- `threedfin_dbh`, `threedfin_hmax` - 3DFIN detection results
- `field_KoealaId` - Plot identifier
- `treesp`, `field_PuuLk` - Tree species codes
- `distance_to_threedfin_cm` - Matching distance threshold
- Various volume, basal area, and point cloud metrics

#### 3. `threedfin_results_A.m`
- **Size:** ~40 KB (40,213 bytes)
- **Type:** MATLAB script (929 lines)
- **Purpose:** Comprehensive tree measurement analysis and validation

## Technical Details

### Data Processing Workflow

The `threedfin_results_A.m` script follows this workflow:

1. **Data Preprocessing**
   - Loads reference field data from external CSV files
   - Filters trees based on plot IDs (44 specific plots)
   - Removes fallen and clipped trees from analysis
   - Filters trees < 5cm DBH
   - Excludes tree class 8 (PuuLk == 8)
   - Removes trees marked as outside plots (Ulkona == 1)

2. **Tree Detection Analysis**
   - Loads combined TLS + 3DFIN dataset
   - Converts DBH units (cm)
   - Filters by distance threshold (< 5cm matching distance)
   - Applies logarithmic height/DBH outlier detection using 2-sigma rule

3. **Statistical Analysis** (Multiple Thresholds: 5, 7, 10 cm)
   - **RMSE** (Root Mean Square Error)
   - **Bias** (systematic error)
   - **Correlation** (Pearson's r and R²)
   - **Detection rates** (segmentation, measurement, finding ratios)

4. **Plot-Level Aggregation**
   - Mean, std, min, max calculations per plot
   - Basal area calculations
   - Stems per hectare
   - Volume estimation using Laasasenaho models
   - Dominant height (Hdom) calculations
   - Basal area weighted means (Dg, Hg)

5. **Output Generation**
   - Plot-level summary statistics
   - Tree-level correlation tables
   - Residual calculations

### Tree Species Codes

The script uses Finnish forestry tree species codes:

| Code | Species (Finnish) | Species (English) |
|------|-------------------|-------------------|
| 1    | Mänty             | Scots Pine        |
| 2    | Kuusi             | Norway Spruce     |
| 3    | Rauduskoivu       | Silver Birch      |
| 4    | Hieskoivu         | Downy Birch       |
| 5    | Haapa             | Aspen             |
| 6    | Harmaaleppä       | Grey Alder        |
| 7    | Tervaleppä        | Common Alder      |
| 8    | Muu havupuu       | Other Conifer     |
| 9    | Muu lehtipuu      | Other Deciduous   |
| 10-29 | Various species  | See script lines 691-703 |

### Volume Calculation Models

**Laasasenaho Volume Models** (lines 706-740):

- **Pine (Mänty):** Different coefficients for DBH and height
- **Spruce (Kuusi):** Species-specific allometric equations
- **Birch and others:** Combined deciduous model
- Trees < 1.3m height: Volume = 0
- Output converted from dm³ to m³

### Key Variables and Metrics

**Tree-Level:**
- `DBH` / `field_D` - Diameter at Breast Height (cm)
- `H` / `field_H` - Tree height (m)
- `ba` - Basal area (m²)
- `vol` / `volumeLaasasen` - Tree volume (m³)

**Plot-Level:**
- `Dg` - Basal area weighted mean diameter
- `Hg` - Basal area weighted mean height
- `Hdom` - Dominant height (mean of 100 largest trees/ha)
- `ba_per_ha` - Basal area per hectare (m²/ha)
- `stems_per_ha` - Tree density (trees/ha)
- `vol_per_ha` - Volume per hectare (m³/ha)

**Quality Metrics:**
- `tree_finding_ratio` - Detection rate (%)
- `RMSE` - Root Mean Square Error
- `Bias` - Mean difference (predicted - observed)
- `r` - Pearson correlation coefficient
- `r2` - Coefficient of determination

## Development Workflows

### Working with This Repository

#### Prerequisites
- MATLAB (any recent version supporting table operations)
- Access to external data files referenced in the script:
  - `Lukupuu_maastolaserkohteet.csv`
  - `kuviot_kessu.csv`
  - These are hardcoded paths to a Windows user directory

#### Common Tasks

**1. Running the Analysis**
```matlab
% Open MATLAB and navigate to repository directory
cd /path/to/silver-bassoon

% Run the main script
threedfin_results_A
```

**2. Modifying Input Data**
- Update CSV file path on line 92
- Adjust plot inclusion list (lines 14-20)
- Modify fallen/clipped tree IDs (lines 26-41)

**3. Changing Analysis Parameters**
- Distance threshold: line 96 (`distance_to_threedfin_cm < 5`)
- DBH thresholds: line 207 (`thresholds = [5, 7, 10]`)
- Outlier detection: line 138 (2-sigma rule)

**4. Adding New Trajectory/Walk**
- Modify line 5: `traj = {'G'};` to include new identifiers
- Script loops through trajectory codes for batch processing

### Path Dependencies

**Critical:** The script contains hardcoded Windows paths that must be updated:

```matlab
% Lines to modify for different environments:
Line 24: 'C:/Users/laurliik/OneDrive - ...'
Line 82: 'C:\Users\laurliik\OneDrive - ...'
Line 92: 'C:\Users\laurliik\Riegl Scans\3Dfin_A\...'
```

**For Linux/Mac users:** Replace backslashes with forward slashes and update paths.

## Key Conventions

### Code Style
- Variable naming: `snake_case` for most variables
- Structure fields: `camelCase` for nested structures
- Comments: Mix of English and Finnish
- Line length: No strict limit

### Data Filtering Rules
1. **Spatial filtering:** Trees must be within 5cm of 3DFIN match
2. **Size filtering:** Only trees ≥ 5cm DBH
3. **Class filtering:** Exclude tree class 8
4. **Boundary filtering:** Remove trees marked "Ulkona" (outside)
5. **Outlier filtering:** 2-sigma rule on height-DBH relationship

### Output Structures

The script creates MATLAB structures:
- `whole_population` - Tree-level statistics by threshold
- `plot_levels` - Plot-level summaries by walk
- `stats_plot_levels` - Aggregated plot statistics
- `correlation_tables` - Correlation matrices
- `All_measured_trees` - Complete tree datasets

## AI Assistant Guidelines

### When Analyzing This Code

1. **Understand the Domain Context**
   - This is forestry research analyzing laser scanning accuracy
   - Field measurements are the "ground truth"
   - TLS/3DFIN measurements are being validated

2. **Data Dependencies**
   - The CSV file is essential - don't delete or modify without backup
   - External CSV files are required but not in repository
   - Script expects specific column names - case-sensitive

3. **Be Cautious With:**
   - Changing statistical thresholds (affects research validity)
   - Modifying species codes or volume models (scientifically validated)
   - Altering hardcoded tree IDs (manually verified data)

4. **Safe Modifications:**
   - Adding visualization code
   - Exporting results to different formats
   - Creating summary reports
   - Improving code readability/documentation

### When Helping Users

1. **Common Issues:**
   - Path errors (Windows vs. Unix)
   - Missing external data files
   - MATLAB version compatibility (table functions)
   - Memory issues with large point cloud datasets

2. **Suggested Improvements:**
   - Parameterize hardcoded paths
   - Add input validation
   - Create configuration file for parameters
   - Add progress indicators for long operations
   - Generate visualization outputs

3. **Documentation Gaps:**
   - No README present
   - Limited inline comments
   - Mixed language comments (English/Finnish)
   - No usage examples

## Research Context

### Study Design
- **Study type:** Terrestrial Laser Scanning validation
- **Location:** Finland (based on coordinate system and language)
- **Plots:** 44 forest inventory plots
- **Methods:** Comparing automated tree detection (3DFIN) vs. manual field measurements
- **Coordinate system:** UTM (likely UTM zone 35N based on coordinates)

### Scientific Standards
- Uses established Finnish forest mensuration methods
- Laasasenaho (1982) volume equations - standard in Finland
- Dominant height based on 100 largest trees per hectare
- Basal area weighted means for stand-level estimates

## Version Control

### Git Workflow
- **Current branch:** `claude/claude-md-miinrpq4wxt9pr5j-01VcrCRXbKLtijt9brv1na9J`
- **Main branch:** Not specified (likely `main` or `master`)
- **Branch naming:** Uses `claude/` prefix for AI-generated branches

### Recent History
```
190f9c6 - Merge pull request #1 from LLiiko/LLiiko-patch-1
3618e54 - Add files via upload
095a973 - Create Train
```

### Commit Guidelines
- Keep data files separate from code changes
- Document statistical parameter changes
- Note any modifications to filtering criteria
- Reference plot IDs when updating tree lists

## Future Development

### Potential Enhancements
1. **Modularity:** Break script into functions
2. **Configuration:** External config file for parameters
3. **Visualization:** Plot generation for results
4. **Automation:** Batch processing for multiple datasets
5. **Documentation:** Add README with research context
6. **Validation:** Unit tests for calculations
7. **Internationalization:** Translate Finnish comments
8. **Performance:** Optimize loops and memory usage

### Data Management
- Consider adding `.gitignore` for large data files
- Document expected data format/schema
- Create example/test datasets
- Add data validation checks

## Contact & Attribution

**Repository Owner:** LLiiko
**Platform:** GitHub
**License:** Not specified (add LICENSE file)

---

*This CLAUDE.md was generated automatically by AI analysis. Last updated: 2025-11-28*
*For questions about the research methodology, consult the repository owner.*
