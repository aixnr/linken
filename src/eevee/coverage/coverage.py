"""eevee/coverage/coverage.py
Input: headerless tab-separated (.tsv) output from samtools mpileup containing 4 columns
  see ./src/hermes/process/templates/coverage.sh

Columns:
  1. Chromosome name
  2. Position
  3. Base
  4. Number of reads
  5. Variant frequency
"""


import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter
from typing import List


class Coverage:
    def __init__(self, sample: str, path_tsv: str):
        _column_names = ["Chr", "Position", "Base", "Coverage"]
        _df = pd.read_csv(path_tsv, sep="\t", header=None, names=_column_names)

        self.df = _df
        self.sample = sample

    def plot_individual_coverage(self, color="DarkSlateGray"):
        """
        First, check how many unique entries (first column), then loop 1 plot for each unique entry
        """
        _df = self.df
        unique_chr: List[str] = _df["Chr"].unique()

        for chr in unique_chr:
            _df_chr = _df.query(f"Chr == '{chr}'")
            
            fig, ax = plt.subplots()
            plot_coverage(ax=ax, df=_df_chr, color=color)

            ax.set_title(f"{self.sample} - {chr}", loc="left", fontsize=9)

            FIG_PATH = f"output_plots/{self.sample} - {chr}.jpg"
            fig.savefig(FIG_PATH, dpi=150, format="jpg", bbox_inches="tight")

    def plot_individual_coverage_variant(self, variant_table: str):
        _df = self.df

        _column_map = {
            "CHROM": "Chr",
            "POS": "Position",
            "REF": "Reference",
            "ALT": "Alternate",
            "AF": "Frequency"
        }

        _df_variant = pd.read_csv(variant_table, sep="\t").rename(columns=_column_map)
        _df_merged = _df.merge(_df_variant, how="inner", on=["Chr", "Position"])

        _df = self.df
        unique_chr: List[str] = _df["Chr"].unique()

        for chr in unique_chr:
            _df_chr = _df.query(f"Chr == '{chr}'")
            _merged_chr = _df_merged.query(f"Chr == '{chr}'")

            # fig, ax = plt.subplots()
            fig = plt.figure()
            gs = fig.add_gridspec(ncols=1, nrows=6)
            ax1 = fig.add_subplot(gs[0:4, :]) # plt.Axes for coverage
            ax2 = fig.add_subplot(gs[4:6, :])   # plt.Axes for allelic frequency

            plt.subplots_adjust(hspace=0.4)

            plot_coverage(ax=ax1, df=_df_chr)
            plot_variant_bar(ax=ax1, df=_merged_chr)
            plot_allelic_frequency(ax=ax2, df=_merged_chr)

            ax1.set_xticklabels([])
            ax1.set_title(f"{self.sample} - {chr}", loc="left", fontsize=9)

            variant_cutoff = 0.01
            ax2.fill_between(x=_df_chr["Position"], y1=variant_cutoff, y2=0, color="silver", alpha=0.5)
            ax2.axhline(y=variant_cutoff, linewidth=0.5, color="silver")

            max_x = _df_chr["Position"].max() + 20
            ax1.set_xlim(0, max_x)
            ax2.set_xlim(0, max_x)

            FIG_PATH = f"output_plots/{self.sample} - {chr} var.jpg"
            fig.savefig(FIG_PATH, dpi=150, format="jpg", bbox_inches="tight")


def plot_coverage(ax: plt.Axes, df: pd.DataFrame, color="DarkSlateGray"):
    max_y = df["Coverage"].max() * 1.25 + 25

    ax.plot(df["Position"], df["Coverage"], color=color)
    ax.fill_between(df["Position"], y1=df["Coverage"], y2=0, color=color, alpha=0.5)

    ax.set_ylim(0, max_y)
    ax.spines[["top", "right"]].set_visible(False)
    ax.set_ylabel("n of reads")


def plot_allelic_frequency(ax: plt.Axes, df: pd.DataFrame, color="Tomato"):
    for _, row in df.iterrows():
        ax.vlines(x=row["Position"], ymin=0, ymax=row["Frequency"], color=color)
        
    ax.set_yscale("log", base=10)
    ax.yaxis.set_major_formatter(ScalarFormatter(useOffset=False)) 
    ax.ticklabel_format(style="plain", axis="y") 

    ax.set_ylim(0.001, 1.316)
    ax.set_yticks([0.001, 0.01, 0.1, 1], [0.001, 0.01, 0.1, 1], fontsize=7)
    ax.set_ylabel("frequency")
    ax.set_xlabel("nucleotide position")


def plot_variant_bar(ax: plt.Axes, df: pd.DataFrame, color="Tomato"):
    for _, row in df.iterrows():
        ax.vlines(x=row["Position"], ymin=0, ymax=row["Coverage"] * 1.1, color=color)

        nuc_ref, nuc_pos, nuc_alt = row['Reference'], row['Position'], row['Alternate']
        nuc_y_height = row["Coverage"] * 1.15
        nuc_marker = f"{nuc_ref}{nuc_pos}{nuc_alt}"

        ax.text(s=nuc_marker, x=nuc_pos, y=nuc_y_height, ha="center",
                color=color, fontsize=7, rotation=90)
