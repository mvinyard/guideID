import pandas as pd

def _annotate_gene_feature(df, feature):

    """"""

    features = df.loc[df["feature"] == feature][["start", "end"]].values
    return pd.DataFrame(data=features, columns=["gene_feature.start", "gene_feature.end"])