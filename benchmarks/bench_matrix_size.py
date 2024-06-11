"""
Comparison of storage size of Phylo2Mat vs. Newick strings
"""

import sys

from argparse import ArgumentParser

import numpy as np
import pandas as pd

from tqdm import tqdm

from benchmarks.plot import plot_sizes
from phylo2vec.matrix import to_newick
from phylo2vec.utils import sample

MIN_LEAVES = 5
MAX_LEAVES = 10000
STEP_LEAVES = 50

tests = [
    "Phylo2Mat (structured)",
    "Phylo2Mat (float16)",
    "Phylo2Mat (object)",
    "Phylo2Mat (string)",
    "Newick (string with BLs rounded to 6 digits)",
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
    sizes = {test: np.zeros((len(all_leaves),), dtype=np.int32) for test in tests}

    for i, n_max in tqdm(enumerate(all_leaves), total=len(all_leaves)):
        v = sample(n_max)
        bls = np.random.uniform(low=0, high=1, size=(n_max - 1, 2)).astype(np.float32)

        m = np.concatenate([v[:, None], bls.astype(np.float16)], axis=1)

        m_obj = np.concatenate(
            [v[:, None], bls.astype(np.float16)], axis=1, dtype=object
        )

        m_str = pd.DataFrame(m_obj).to_csv(index=False)

        m_structured = np.array(
            [(v[i], bls[i, 0], bls[i, 1]) for i in range(n_max - 1)],
            dtype=[("v", np.int16), ("bl1", np.float16), ("bl2", np.float16)],
        )

        newick_six = to_newick(v, bls.round(6).astype(str))

        newick = to_newick(v, bls.astype(str))

        sizes["Phylo2Mat (float16)"][i] = sys.getsizeof(m)
        sizes["Phylo2Mat (object)"][i] = sys.getsizeof(m_obj)
        sizes["Phylo2Mat (structured)"][i] = sys.getsizeof(m_structured)
        sizes["Phylo2Mat (string)"][i] = sys.getsizeof(m_str)
        sizes["Newick (string with BLs rounded to 6 digits)"][i] = sys.getsizeof(
            newick_six
        )
        sizes["Newick (string)"][i] = sys.getsizeof(newick)

    sizes_df = pd.DataFrame(sizes)
    sizes_df["n_leaves"] = all_leaves
    sizes_df = sizes_df.melt(
        id_vars="n_leaves", var_name="format", value_name="size [B]"
    )
    sizes_df["size [kB]"] = sizes_df["size [B]"].div(1000)

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
