
import pandas as pd

from ._number_features import _number_features


def _format_feature_df(feature_df, feature, reverse_strand):
    
    """
    Format the feature dataframe such that the columns are compatible with that of the guideID search module. 
    
    Parameters:
    -----------
    feature_df
    
    feature
    
    Returns:
    --------
    df
        formatted DataFrame
    
    Notes:
    ------
    """

    df_ = feature_df.copy()
    start = df_.filter(regex="tart")
    end = df_.filter(regex="nd")
    df = pd.concat([df_['Chromosome'], start, end], axis=1)
    df.columns = ["Chromosome", "gene_feature.start", "gene_feature.end"]

    df["gene_feature"] = feature
    df["length"] = df["gene_feature.end"] - df["gene_feature.start"]
    
    df = _number_features(df, feature, reverse_strand)
    
    return df