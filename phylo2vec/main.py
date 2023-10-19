from argparse import ArgumentParser

import numpy as np

from phylo2vec.to_newick import to_newick
from phylo2vec.utils import check_v


def parse_args():
    """Parse main arguments."""
    parser = ArgumentParser(prog="Phylo2Vec", description="Phylo2vec arguments")

    parser.add_argument(
        "--toNewick",
        default=False,
        action="store_true",
        help="(bool) Convert to Newick format. Example input: 0,1,4",
    )

    parser.add_argument(
        "--toVector",
        default=False,
        action="store_true",
        help='(bool) Convert to integer vector. Example input: "(((2,1)4,0)5,3)6;"',
    )

    parser.add_argument(
        "--withMapping",
        default=False,
        action="store_true",
        help="(bool) For Newicks that do not only contain digits, to use with toVector.\n"
        'Example input: "(((((((tip_0:1.44,tip_1:1.44)8042:0.46,(tip_2:1.5,tip_3:1.5)8043:0.4)8044:0.3,(tip_4:1.51,tip_5:1.51)8045:0.69)8046:0.4,tip_6:2.6)8047:1.05,tip_7:3.65)8048:0.5,(((tip_8:0.72,tip_9:0.72)8049:0.28,tip_10:1)8050:1.56,tip_11:2.56)8051:1.59)8052:1.96,tip_12:6.11)8053:0;"',
    )

    parser.add_argument(
        "--num_leaves",
        type=int,
        help="(int) Number of leaves (optional, but recommended when using toVector)",
    )

    parser.add_argument(
        "--verbose",
        default=False,
        action="store_true",
        help="(bool) If True, display intermediate results",
    )

    parser.add_argument(
        "--input",
        required=True,
        type=str,
        help="(str, required) Input to toNewick or toVector",
    )

    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()

    if args.toNewick:
        v = np.fromiter(map(int, args.input.split(",")), dtype=np.int64)

        check_v(v)

        print(to_newick(v))

    elif args.toVector:
        pass
    else:
        raise TypeError("Expected --toNewick or --toVector.")
