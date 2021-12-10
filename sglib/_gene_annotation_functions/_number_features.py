import numpy as np
import pandas as pd


def _number_features(guide_df, feature="exon", reverse_strand=False):

    """"""

    exon_starts = guide_df["gene_feature.start"].unique()

    if not reverse_strand:
        exons = np.arange(len(exon_starts))
    else:
        exons = np.flip(np.arange(len(exon_starts)))

    df_ = pd.DataFrame(
        data=np.stack([exon_starts, exons]).T, columns=["gene_feature.start", feature]
    )

    guide_df = pd.merge(guide_df, df_, on="gene_feature.start")

    return guide_df