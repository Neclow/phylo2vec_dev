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


def remove_annotations_from_newick(newick):
    """Remove branch lengths annotations from a Newick string

    Example: "(((2:0.02,1:0.01),0:0.041),3:1.42);" --> "(((2,1),0),3);"

    Parameters
    ----------
    newick : str
        Newick representation of a tree
    """
    return re.sub(ANNOTATION_PATTERN, "", newick)
