import vinplots
import matplotlib.pyplot as plt


def _calculate_library_summary_metrics(guide_df, feature="feature"):

    """
    For a CDS, this would be "exon". 
    """
    
    guide_df["length"] = abs(guide_df["Start"] - guide_df["End"])
    groupby_region = guide_df.groupby(feature)
    region_lengths = groupby_region.describe()["length"]["mean"]
    guides_per_region = groupby_region.count()["protospacer.start"].values

    guide_per_bp = guides_per_region / region_lengths
    bp_per_guide = 1 / guide_per_bp

    regions = region_lengths.index

    return regions, guides_per_region, guide_per_bp, bp_per_guide


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


def _plot_guide_distribution(guide_df, feature):

    """"""

    fig, axes = _construct_plot()
    titles = ["n_guides / {}".format(feature), "n_guides / bp / {}".format(feature), "n_bp / guide"]

    (
        regions,
        guides_per_region,
        guide_per_bp,
        bp_per_guide,
    ) = _calculate_library_summary_metrics(guide_df, feature)
    
    plot_guide_stats = [guides_per_region, guide_per_bp, bp_per_guide]

    for i, ax in enumerate(axes):
        ax.set_title(titles[i])
        ax.grid(True, zorder=0, alpha=0.2)
        ax.bar(x=regions, height=plot_guide_stats[i], color="k", zorder=3)