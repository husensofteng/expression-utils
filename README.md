# Log2FC Contrast Pair Correlation Plot

A tool for visualizing log2 fold-change correlations between differential expression contrasts.

## Installation

```bash
pip install pandas numpy matplotlib scipy
```

## Usage

```bash
python plot_log2fc_contrast_pair_correlation.py [OPTIONS]
```

### Parameters

| Parameter | Description | Default |
|-----------|-------------|---------|
| `--deg-file` | DEG file in `Label=path` format. Can be specified multiple times. | Built-in defaults |
| `--deg-dir` | Directory containing DEG files | `deg_files` |
| `--output-dir` | Output directory for plots | `.` (current directory) |

### Input File Format

Tab-delimited files with the following columns:
- `gene_id` - Gene identifier
- `log2FoldChange` - Log2 fold change value
- `padj` - Adjusted p-value

### Examples

```bash
# Using default files
python plot_log2fc_contrast_pair_correlation.py

# Specify custom DEG files
python plot_log2fc_contrast_pair_correlation.py \
    --deg-file "Contrast1=treatment1_vs_control.tsv" \
    --deg-file "Contrast2=treatment2_vs_control.tsv"

# Specify directories
python plot_log2fc_contrast_pair_correlation.py \
    --deg-dir /path/to/deg_files \
    --output-dir /path/to/output
```

## Output

- `log2fc_comparison_plot.svg` - Vector format plot
- `log2fc_comparison_plot.pdf` - PDF format plot (300 DPI)

## Plot Features

- **Color coding by significance** (padj < 0.05):
  - Green: Significant in both contrasts
  - Orange: Significant in x-axis contrast only
  - Blue: Significant in y-axis contrast only
  - Gray: Not significant in either contrast

- **Point size**: Scaled by mean adjusted p-value (-log10)

- **Correlation lines**: Pearson R shown for "Both": genes that are significannt in both contrasts and "Neither" for genes that are not significant in either contrast.

- **Outlier handling**: X-axis values < -10 are displayed in a compressed gap region
