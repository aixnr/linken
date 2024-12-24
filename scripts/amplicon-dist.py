#!/pyvenv/bin/python

"""
Output a histogram of amplicon length.
Was designed for ONT sequence reads, but can be used for other platforms.

Usage:
  python amplicon-dist[.py] --gz filename.fastq.gz
  python amplicon-dist[.py] --fq filename.fastq
"""

import gzip
import argparse
import statistics
import matplotlib.pyplot as plt
from pathlib import Path


def check_if_exists(location: str) -> None:
    if not Path(location).exists():
        print("")
        print(f"  File '{location}' is absent, aborting...")
        exit(1)


def read_fastq(handle) -> list:
    lines = handle.readlines()
    file_lines_total = len(lines)

    line_skip_for_basecall = 4
    skip_range = range(1, file_lines_total, line_skip_for_basecall)
    skip_range_list = [x for x in skip_range]

    bases_list_length = [len(lines[i].decode().strip()) for i in skip_range_list]
    bases_list_length_200bp = [x for x in bases_list_length if x >= 200]

    # 200bp = to exclude primer dimer
    return bases_list_length_200bp


def generate_histogram_length(archive: str = None, fastq: str = None) -> None:
    if archive:
        _name = archive
        with gzip.open(archive, "rb") as f:
            _amplicon_length = read_fastq(handle=f)
    elif fastq:
        _name = fastq
        with open(fastq, "rb") as f:
            _amplicon_length = read_fastq(handle=f)

    _mode_count = statistics.mode(_amplicon_length)

    fig, ax = plt.subplots()
    ax.hist(_amplicon_length, bins="auto", alpha=0.7, color="DodgerBlue", edgecolor="DarkSlateGray")

    _max_y = ax.get_ylim()[1]
    for _y in range(0, int(_max_y), 50):
        ax.axhline(y=_y, color="silver", linewidth=0.5, alpha=0.25, zorder=-100)

    ax.set_xlabel("Distribution of amplicon length [bp]")
    ax.set_ylabel("Count [n]")
    ax.set_title(f"{_name} (mode: {_mode_count:,} bp)", loc="left", fontsize=10)
    ax.spines[["top", "right"]].set_visible(False)

    fig.savefig(f"{_name}.jpg", format="jpg", dpi=150, bbox_inches="tight")


# ------------------------------------------------------------------------------
# Argparse CLI
def cli() -> None:
    parser = argparse.ArgumentParser()

    parser.add_argument("--gz", "-g", required=False, type=str, help="Specify the location of the '.fastq.gz' archive file.")
    parser.add_argument("--fq", "-f", required=False, type=str, help="Specify the location of the '.fastq' file.")

    args = parser.parse_args()

    if args.gz:
        check_if_exists(args.gz)
        generate_histogram_length(archive=args.gz)
    elif args.fq:
        check_if_exists(args.fq)
        generate_histogram_length(fastq=args.fq)
    else:
        parser.print_help()


# ------------------------------------------------------------------------------
# main
if __name__ == "__main__":
    cli()
