import pandas as pd

from ._fetch_chromosome_name import _fetch_chromosome_name


def _annotate_gene_feature(df, feature):

    """"""

    features = df.loc[df["feature"] == feature][["start", "end"]].values
    feature_df = pd.DataFrame(data=features, columns=["gene_feature.start", "gene_feature.end"])
    feature_df['Chromosome'] = _fetch_chromosome_name(df)
    
    return feature_df