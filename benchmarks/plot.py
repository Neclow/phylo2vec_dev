"""
Plotting utility functions for the benchmarks
"""

# pylint: disable=protected-access

import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns


def clear_axes(
    ax=None, top=True, right=True, left=False, bottom=False, minorticks_off=True
):
    """A more forcing version of sns.despine.

    Parameters
    ----------
    ax : matplotlib axes, optional
        Specific axes object to despine.
    top, right, left, bottom : boolean, optional
        If True, remove that spine.
    minorticks_off: boolean, optional
        If True, remove all minor ticks
    """
    if ax is None:
        axes = plt.gcf().axes
    else:
        axes = [ax]

    for ax_i in axes:
        sns.despine(ax=ax_i, top=top, bottom=bottom, left=left, right=right)
        if minorticks_off:
            ax_i.minorticks_off()

        ax_i.tick_params(
            top=not top,
            bottom=not bottom,
            left=not left,
            right=not right,
            labelleft=not left,
            labelright=not right,
            labeltop=not top,
            labelbottom=not bottom,
        )


def set_size(width, layout="h", fraction=1):
    """Set figure dimensions in inches to avoid scaling in LaTeX.

    Adapted from: https://jwalton.info/Embed-Publication-Matplotlib-Latex/

    Parameters
    ----------
    width: float
        Document textwidth or columnwidth in pts
        Report: 390 pt
    layout: string
        h: horizontal layout
        v: vertical layout
        s: square layout
    fraction: float, optional
        Fraction of the width which you wish the figure to occupy

    Returns
    -------
    fig_dim: tuple
            Dimensions of figure in inches
    """
    # Width of figure (in pts)
    fig_width_pt = width * fraction

    # Convert from pt to inches
    inches_per_pt = 1 / 72.27

    # Golden ratio for aesthetic figures
    # https://disq.us/p/2940ij3
    golden_ratio = (5**0.5 - 1) / 2

    # Figure width in inches
    fig_width_in = fig_width_pt * inches_per_pt
    # Figure height in inches
    if layout == "h":
        fig_height_in = fig_width_in * golden_ratio
    elif layout == "v":
        fig_height_in = fig_width_in / golden_ratio
    elif layout == "s":
        fig_height_in = fig_width_in

    fig_dim = (fig_width_in, fig_height_in)

    return fig_dim


def plot_sizes(tests, sizes_df, output_pdf, no_latex, show_plot):
    if not no_latex:
        plt.rcParams.update(
            {
                "font.serif": ["sffamily"],
                "figure.dpi": "100",
                "font.size": 9,
                "text.usetex": True,
            }
        )

        matplotlib.rc("text.latex", preamble=r"\usepackage{amsmath}")

    _, ax = plt.subplots(1, 1, figsize=set_size(290, "h"))

    palette = dict(zip(tests, [*sns.color_palette("OrRd_r", n_colors=3).as_hex(), "k"]))

    sns.scatterplot(
        x="n_leaves",
        y="Size [kB]",
        hue="format",
        data=sizes_df,
        s=4,
        palette=palette,
        ax=ax,
    )

    ax.legend(title="")
    ax.xaxis.get_major_formatter()._usetex = False
    ax.yaxis.get_major_formatter()._usetex = False

    clear_axes()

    plt.savefig(output_pdf, bbox_inches="tight", pad_inches=0)

    print(f"Plot saved at {output_pdf}")

    if show_plot:
        plt.show()
