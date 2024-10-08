"""
Comparison of storage size of Phylo2Vec vectors vs. Newick strings
"""

import sys

from argparse import ArgumentParser

import numpy as np
import pandas as pd

from tqdm import tqdm

from benchmarks.plot import plot_sizes
from phylo2vec.base import to_newick
from phylo2vec.utils import sample

MIN_LEAVES = 5
MAX_LEAVES = 10000
STEP_LEAVES = 50

tests = [
    "Phylo2Vec (int16)",
    "Phylo2Vec (int32)",
    "Phylo2Vec (string)",
    "Newick (string)",
]


def parse_args():
    """Parse optional arguments."""
    parser = ArgumentParser(description="Tree format size benchmark tool")
    parser.add_argument(
        "--show_plot",
        action="store_true",
        help="Show plot output with matplotlib.",
    )
    parser.add_argument(
        "--no-latex",
        action="store_false",
        help="Do not use LaTeX fonts for math in matplotlib",
    )
    parser.add_argument(
        "--output-file",
        type=str,
        default="bench_size",
        help="Output file name",
        required=True,
    )

    return parser.parse_args()


def compute_sizes(all_leaves, output_csv):
    # Pre-allocate size arrays for each test
    sizes = {test: np.zeros((len(all_leaves),), dtype=np.int32) for test in tests}

    # Compute sizes
    for i, n_leaves in tqdm(enumerate(all_leaves), total=len(all_leaves)):
        v = sample(n_leaves)
        newick = to_newick(v)

        sizes["Phylo2Vec (int16)"][i] = sys.getsizeof(v)
        sizes["Phylo2Vec (int32)"][i] = sys.getsizeof(v.astype(np.int32))
        sizes["Phylo2Vec (string)"][i] = sys.getsizeof(",".join(map(str, v)))
        sizes["Newick (string)"][i] = sys.getsizeof(newick)

    # Make a DataFrame and convert to long format
    sizes_df = pd.DataFrame(sizes)
    sizes_df["n_leaves"] = all_leaves

    sizes_df = sizes_df.melt(
        id_vars="n_leaves", var_name="format", value_name="Size [B]"
    )

    sizes_df["Size [kB]"] = sizes_df["Size [B]"].div(1000)

    sizes_df.to_csv(output_csv, index=False)

    return sizes_df


def main():
    """Main script"""
    args = parse_args()

    output_csv = f"benchmarks/res/{args.output_file}.csv"
    output_pdf = f"benchmarks/img/{args.output_file}.pdf"

    all_leaves = np.arange(MIN_LEAVES, MAX_LEAVES, STEP_LEAVES)

    sizes_df = compute_sizes(all_leaves, output_csv)

    print(f"Data saved at {output_csv}")

    # Make the plot
    plot_sizes(
        tests, sizes_df, output_pdf, no_latex=args.no_latex, show_plot=args.show_plot
    )


if __name__ == "__main__":
    main()
