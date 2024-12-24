#!/pyvenv/bin/python

"""
A simple tool for wrapping text to a predefined width.
Default is 60 characters long.
"""

import textwrap as tr
import argparse


def wrapping(length: int) -> None:
    seq = input("Paste the sequence here: ")
    _seq = seq.strip().replace(" ", "")

    _tr = tr.fill(_seq, width=length)
    print("")
    print(_tr)
    print("")
    print(f"Length: {len(_seq)}\n")


def cli():
    parser = argparse.ArgumentParser()
    parser.add_argument("--length", "-l", type=int, required=False, default=60, help="Length per line.")

    args = parser.parse_args()
    wrapping(length=args.length)


if __name__ == "__main__":
    cli()
