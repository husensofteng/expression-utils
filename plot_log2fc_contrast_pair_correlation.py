#!/usr/bin/env python3
"""
Create log2FC-vs-log2FC comparison plots for DE genes across conditions.
Defaults are provided, but inputs can be overridden via CLI arguments.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
from itertools import combinations
import os
import argparse
from functools import reduce

DEFAULT_DEG_DIR = "deg_files"
DEFAULT_OUTPUT_DIR = "."
DEFAULT_DEG_FILES = {
    "contrast1": os.path.join(DEFAULT_DEG_DIR, "contrast1.tsv"),
    "contrast2": os.path.join(DEFAULT_DEG_DIR, "contrast2.tsv"),
    "contrast3": os.path.join(DEFAULT_DEG_DIR, "contrast3.tsv"),
}

# Padj threshold for significance
PADJ_THRESHOLD = 0.05


def load_deg_data(filepath, label):
    """Load DEG data and extract relevant columns."""
    df = pd.read_csv(str(filepath), sep="\t")
    # Standardize column names
    df = df.rename(columns={
        "gene_id": "gene_id",
        "log2FoldChange": "log2FC",
        "padj": "padj",
        "pvalue": "pvalue"
    })
    # Select only necessary columns and rename log2FC with label
    result = df[["gene_id", "log2FC", "padj"]].copy()
    result = result.rename(columns={
        "log2FC": f"log2FC_{label}",
        "padj": f"padj_{label}"
    })
    return result


def significance_category(padj_x, padj_y, x_label, y_label):
    """Classify gene based on significance in each dataset."""
    sig_x = (not pd.isna(padj_x)) and (padj_x < PADJ_THRESHOLD)
    sig_y = (not pd.isna(padj_y)) and (padj_y < PADJ_THRESHOLD)

    if sig_x and sig_y:
        return "Both"
    elif sig_x and not sig_y:
        return f"{x_label} only"
    elif not sig_x and sig_y:
        return f"{y_label} only"
    else:
        return "Neither"


def get_point_size(padj, min_size=10, max_size=150):
    """Scale point size with significance (-log10 for size only)."""
    if pd.isna(padj) or padj <= 0:
        return min_size
    neg_log10 = -np.log10(padj)
    size = min_size + (neg_log10 / 10) * (max_size - min_size)
    return min(max(size, min_size), max_size)


def create_comparison_plot(data, x_col, y_col, x_label, y_label, ax):
    """Create a single log2FC comparison scatter plot."""
    plot_data = data[[x_col, y_col, f"padj_{x_label}", f"padj_{y_label}"]].dropna()

    plot_data["mean_padj"] = (
        plot_data[f"padj_{x_label}"] + plot_data[f"padj_{y_label}"]
    ) / 2
    # Categorize by significance in each dataset
    plot_data["sig_cat"] = plot_data.apply(
        lambda row: significance_category(
            row[f"padj_{x_label}"], row[f"padj_{y_label}"], x_label, y_label
        ), axis=1
    )
    plot_data["size"] = plot_data["mean_padj"].apply(get_point_size)

    # Four colors based on significance pattern
    colors = {
        "Both": "#2E8B57",              # Green - significant in both
        "Neither": "#808080",            # Gray - not significant in either
        f"{x_label} only": "#E69F00",    # Orange - significant in x only
        f"{y_label} only": "#56B4E9",    # Sky blue - significant in y only
    }

    category_order = [
        "Neither",
        f"{x_label} only",
        f"{y_label} only",
        "Both",
    ]

    # Count genes per category for legend
    category_counts = {cat: len(plot_data[plot_data["sig_cat"] == cat]) for cat in category_order}

    # Check for outliers on x-axis (log2FC < -10)
    outlier_threshold = -10
    has_outliers = (plot_data[x_col] < outlier_threshold).any()

    if has_outliers:
        # Separate outliers from main data
        outlier_data = plot_data[plot_data[x_col] < outlier_threshold]
        main_data = plot_data[plot_data[x_col] >= outlier_threshold]

        # Get x-axis range for main data
        x_min_main = main_data[x_col].min() if len(main_data) > 0 else outlier_threshold
        x_max = plot_data[x_col].max()

        # Add padding
        x_range = x_max - x_min_main
        x_min_main_padded = x_min_main - 0.05 * x_range
        x_max_padded = x_max + 0.05 * x_range

        # Plot outliers at a compressed position (gap region)
        gap_position = x_min_main_padded - 0.08 * x_range

        for cat in category_order:
            subset = outlier_data[outlier_data["sig_cat"] == cat]
            if len(subset) == 0:
                continue
            # Plot outliers at the gap position
            ax.scatter(
                [gap_position] * len(subset),
                subset[y_col],
                s=subset["size"],
                c=colors.get(cat, "#808080"),
                alpha=0.6,
                edgecolors="none",
                marker='<',  # Left-pointing triangle to indicate outliers
            )

        # Plot main data
        for cat in category_order:
            subset = main_data[main_data["sig_cat"] == cat]
            if len(subset) == 0:
                continue
            ax.scatter(
                subset[x_col],
                subset[y_col],
                s=subset["size"],
                c=colors.get(cat, "#808080"),
                alpha=0.6,
                label=cat,
                edgecolors="none",
            )

        # Set x-axis limits with gap
        ax.set_xlim(gap_position - 0.5, x_max_padded)

        # Draw break lines to indicate gap
        break_x = x_min_main_padded - 0.04 * x_range
        y_lim = ax.get_ylim()
        y_range = y_lim[1] - y_lim[0]

        # Draw diagonal break lines
        break_width = 0.02 * x_range
        break_height = 0.02 * y_range
        for y_pos in [y_lim[0] + 0.02 * y_range, y_lim[0] + 0.05 * y_range]:
            ax.plot([break_x - break_width, break_x + break_width],
                   [y_pos - break_height, y_pos + break_height],
                   color='black', linewidth=1, clip_on=False)

        # Add annotation for outliers
        n_outliers = len(outlier_data)
        ax.annotate(f'n={n_outliers}\n(<{outlier_threshold})',
                   xy=(gap_position, y_lim[0] + 0.1 * y_range),
                   fontsize=7, ha='center', va='bottom')

        # Use full data for correlation calculations
        plot_data_for_corr = plot_data
    else:
        # No outliers - plot normally
        for cat in category_order:
            subset = plot_data[plot_data["sig_cat"] == cat]
            if len(subset) == 0:
                continue
            ax.scatter(
                subset[x_col],
                subset[y_col],
                s=subset["size"],
                c=colors.get(cat, "#808080"),
                alpha=0.6,
                label=cat,
                edgecolors="none",
            )
        plot_data_for_corr = plot_data

    ax.axhline(y=0, color="black", linestyle="-", linewidth=0.5)
    ax.axvline(x=0, color="black", linestyle="-", linewidth=0.5)

    # Show correlation only for "Both" and "Neither" categories
    corr_categories = ["Both", "Neither"]
    y_text_positions = [0.95, 0.88]
    for idx, cat in enumerate(corr_categories):
        subset = plot_data_for_corr[plot_data_for_corr["sig_cat"] == cat]
        if len(subset) > 5:
            r, p = stats.pearsonr(subset[x_col], subset[y_col])
            n_genes = len(subset)

            slope, intercept = np.polyfit(subset[x_col], subset[y_col], 1)
            x_line = np.linspace(subset[x_col].min(), subset[x_col].max(), 100)
            y_line = slope * x_line + intercept

            ax.plot(x_line, y_line, color=colors.get(cat, "#808080"),
                   linewidth=1.5, alpha=0.8)

            p_text = "p < 2.2e-16" if p < 2.2e-16 else f"p = {p:.2g}"
            ax.text(
                0.02, y_text_positions[idx],
                f"R = {r:.2f}, {p_text}, n = {n_genes}",
                transform=ax.transAxes,
                fontsize=9,
                color=colors.get(cat, "#808080")
            )

    ax.set_xlabel(f"RNA-seq log2FC\n{x_label}", fontsize=11)
    ax.set_ylabel(f"RNA-seq log2FC\n{y_label}", fontsize=11)

    from matplotlib.lines import Line2D
    legend_elements = [
        Line2D([0], [0], marker='o', color='w', markerfacecolor=colors["Both"],
               markersize=8, label=f"Both (n={category_counts['Both']})", alpha=0.6),
        Line2D([0], [0], marker='o', color='w', markerfacecolor=colors[f"{x_label} only"],
               markersize=8, label=f"{x_label} only (n={category_counts[f'{x_label} only']})", alpha=0.6),
        Line2D([0], [0], marker='o', color='w', markerfacecolor=colors[f"{y_label} only"],
               markersize=8, label=f"{y_label} only (n={category_counts[f'{y_label} only']})", alpha=0.6),
        Line2D([0], [0], marker='o', color='w', markerfacecolor=colors["Neither"],
               markersize=8, label=f"Neither (n={category_counts['Neither']})", alpha=0.6),
    ]

    legend1 = ax.legend(
        handles=legend_elements,
        title=f"Significant\n(padj < {PADJ_THRESHOLD})",
        loc="upper right",
        fontsize=8,
        title_fontsize=8,
        framealpha=0.9
    )
    ax.add_artist(legend1)

    size_legend_elements = []
    for padj_val, label in [(1e-10, "1e-10"), (1e-4, "1e-04"), (1e-2, "1e-02"), (1e-1, "1e-01")]:
        size = get_point_size(padj_val)
        size_legend_elements.append(
            Line2D([0], [0], marker='o', color='w', markerfacecolor='gray',
                   markersize=np.sqrt(size), label=label, alpha=0.6)
        )

    _ = ax.legend(
        handles=size_legend_elements,
        title="Adj. p-value\n(mean)",
        loc="lower right",
        fontsize=7,
        title_fontsize=8,
        framealpha=0.9
    )

    return ax


def parse_args():
    parser = argparse.ArgumentParser(
        description="Create log2FC comparison plots across conditions."
    )
    parser.add_argument(
        "--deg-file",
        action="append",
        help="Label=path for a DEG file (tab-delimited with gene_id, log2FoldChange, padj). "
             "Can be given multiple times. Defaults to built-in set.",
    )
    parser.add_argument(
        "--deg-dir",
        default=DEFAULT_DEG_DIR,
        help=f"Directory containing DEG files (default: {DEFAULT_DEG_DIR}).",
    )
    parser.add_argument(
        "--output-dir",
        default=DEFAULT_OUTPUT_DIR,
        help=f"Output directory for plots (default: {DEFAULT_OUTPUT_DIR}).",
    )
    return parser.parse_args()


def resolve_deg_files(args):
    if not args.deg_file:
        return DEFAULT_DEG_FILES

    deg_files = {}
    for entry in args.deg_file:
        if "=" not in entry:
            raise ValueError(f"Invalid --deg-file format: {entry}. Use label=path.")
        label, path = entry.split("=", 1)
        deg_files[label] = path
    return deg_files


def merge_on_gene_id(dataframes):
    return reduce(lambda left, right: left.merge(right, on="gene_id", how="outer"), dataframes)


def main():
    args = parse_args()
    deg_files = resolve_deg_files(args)

    print("Loading DEG files...")
    deg_data = {
        label: load_deg_data(
            filepath if os.path.isabs(filepath) else os.path.join(args.deg_dir, filepath),
            label
        )
        for label, filepath in deg_files.items()
    }

    print("Merging data...")
    merged = merge_on_gene_id(list(deg_data.values()))
    print(f"Total genes after merge: {len(merged)}")

    labels = list(deg_files.keys())
    pairs = list(combinations(labels, 2))

    n_pairs = len(pairs)
    fig, axes = plt.subplots(1, n_pairs, figsize=(6*n_pairs, 6))
    if n_pairs == 1:
        axes = [axes]

    for idx, (label1, label2) in enumerate(pairs):
        print(f"Creating plot: {label1} vs {label2}...")
        x_col = f"log2FC_{label1}"
        y_col = f"log2FC_{label2}"
        print(len(merged), len(x_col), len(y_col))
        create_comparison_plot(
            merged,
            x_col, y_col,
            label1, label2,
            axes[idx]
        )

        axes[idx].text(
            -0.1, 1.05,
            chr(65 + idx),
            transform=axes[idx].transAxes,
            fontsize=16,
            fontweight='bold'
        )

    plt.tight_layout()

    os.makedirs(args.output_dir, exist_ok=True)
    output_path = os.path.join(args.output_dir, "log2fc_comparison_plot.svg")
    plt.savefig(output_path, bbox_inches="tight", facecolor="white")
    print(f"Saved plot to: {output_path}")

    pdf_path = os.path.join(args.output_dir, "log2fc_comparison_plot.pdf")
    plt.savefig(pdf_path, bbox_inches="tight", facecolor="white", dpi=300)
    print(f"Saved PDF to: {pdf_path}")

    plt.show()


if __name__ == "__main__":
    main()
