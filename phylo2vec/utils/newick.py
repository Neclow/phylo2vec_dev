import re


LEFT_NODE_PATTERN = re.compile(r"\(\b\d+\b")
RIGHT_NODE_PATTERN = re.compile(r",\b\d+\b")
ANNOTATION_PATTERN = re.compile(r":\d+(\.\d+)?")


def find_num_leaves(newick):
    """Calculate the number of leaves in a tree from its Newick

    Parameters
    ----------
    newick : str
        Newick representation of a tree

    Returns
    -------
    int
        Number of leaves
    """
    return len(re.findall(LEFT_NODE_PATTERN, newick)) + len(
        re.findall(RIGHT_NODE_PATTERN, newick)
    )


def create_label_mapping(newick):
    label_mapping = {}

    newick = newick[:-1]  # For ";"

    newick_clean = newick[:-1]  # For ";"

    def do_reduce(newick, newick_clean, j):
        for i, char in enumerate(newick):
            if char == "(":
                open_idx = i + 1
            elif char == ")":
                for child in newick[open_idx:i].split(",", 2):
                    if child not in label_mapping:
                        label_mapping[f"{j}"] = child
                        newick_clean = newick_clean.replace(child, f"{j}")
                        j += 1

                parent = newick[i + 1 :].split(",", 1)[0].split(")", 1)[0]

                newick = newick.replace(
                    newick[open_idx - 1 : i + 1 + len(parent)], f"{j - 1}"
                )

                newick_clean = newick_clean.replace(parent, "")

                return do_reduce(newick, newick_clean, j)

        return newick_clean

    newick_clean = do_reduce(newick, newick_clean, 0)

    return label_mapping, newick_clean + ";"


def remove_annotations(newick):
    """Remove branch lengths annotations from a Newick string

    Example: "(((2:0.02,1:0.01),0:0.041),3:1.42);" --> "(((2,1),0),3);"

    Parameters
    ----------
    newick : str
        Newick representation of a tree
    """
    return re.sub(ANNOTATION_PATTERN, "", newick)


def remove_parent_labels(newick):
    open_idx = -1
    newick_no_parent = newick
    for i, char in enumerate(newick):
        if (char == "," or char == ")" or char == ";") and open_idx != -1:
            newick_no_parent = newick_no_parent.replace(newick[open_idx:i], "")
            open_idx = -1
        if char == ")":
            open_idx = i + 1

    return newick_no_parent
