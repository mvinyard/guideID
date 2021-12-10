def _print_exon_statistics(feature_df):

    groupby_feature = feature_df.groupby("gene_feature")
    return groupby_feature.describe()["length"].T