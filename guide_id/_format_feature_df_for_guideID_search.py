
import pandas as pd


def _format_feature_df_for_guideID_search(feature_df, feature):
    
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
    df = pd.concat([start, end], axis=1)
    df.columns = ["gene_feature.start", "gene_feature.end"]

    df["gene_feature"] = feature
    df["length"] = df["gene_feature.end"] - df["gene_feature.start"]

    return df