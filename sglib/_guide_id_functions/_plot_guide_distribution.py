import vinplots
import matplotlib.pyplot as plt


def _calculate_library_summary_metrics(guide_df):

    """"""
    groupby_exon = guide_df.groupby("exon")
    exon_lengths = groupby_exon.describe()["length"]["mean"]
    guides_per_exon = groupby_exon.count()["protospacer.start"].values

    guide_per_bp = guides_per_exon / exon_lengths
    bp_per_guide = 1 / guide_per_bp

    exons = exon_lengths.index

    return exons, guides_per_exon, guide_per_bp, bp_per_guide


def _construct_plot():

    fig = vinplots.Plot()
    fig.construct(nplots=3, ncols=3, figsize=1.2)

    fig.modify_spines(
        ax="all",
        spines_to_delete=["top", "right"],
        spines_to_move=["bottom", "left"],
        spines_positioning_amount=5,
    )

    axes = [fig.AxesDict[0][0], fig.AxesDict[0][1], fig.AxesDict[0][2]]

    return fig, axes


def _plot_guide_distribution(guide_df):

    """"""

    fig, axes = _construct_plot()
    titles = ["n_guides / exon", "n_guides / bp / exon", "n_bp / guide"]

    (
        exons,
        guides_per_exon,
        guide_per_bp,
        bp_per_guide,
    ) = _calculate_library_summary_metrics(guide_df)
    plot_guide_stats = [guides_per_exon, guide_per_bp, bp_per_guide]

    for i, ax in enumerate(axes):
        ax.set_title(titles[i])
        ax.grid(True, zorder=0, alpha=0.2)
        ax.bar(x=exons, height=plot_guide_stats[i], color="k", zorder=3)