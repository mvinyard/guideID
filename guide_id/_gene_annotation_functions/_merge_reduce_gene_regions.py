import seq_toolkit

def _merge_reduce_gene_regions(gene_df, feature_df, out_path=False, return_df=False):

    """
    Often, when isolating a feature from a GTF for a given gene (ex: exons), there are many overlapping
    isoforms that can be collapsed to avoid redundancy in identifying sgRNAs. This function uses the more
    general `GenomicFeatures.merge()` class from seq_toolkit to accomplish this feature merge-reducing.
    Optionally, this function saves the resulting merged feature set as a tab-delimited .bed file.
    
    Parameters:
    -----------
    gene_df
        subset of the gtf DataFrame obtained using the gtfparse library.
        type: pandas.DataFrame
        
    feature_df
        non-merged-reduced DataFrame of features of interest (subset of the gene_df).
        type: pandas.DataFrame
        
    out_path
        Indicates if merged_feature_df should be saved to disk as a bed file. 
        default: False
        type: bool
        
    return_df:
        Indicates if merged_feature_df should be returned. 
        default: False
        type: bool
    
    Returns:
    --------
    merged_feature_df
        merged-reduced (elimination of overlapping features) DataFrame
        type: pandas.DataFrame
    
    Notes:
    ------
    (1) It's necessary to pass both the gene_df and feature_df in order to get the chromosome name.
    """

    chromosome = _get_chromosome_name(gene_df)
    df_ = feature_df.copy()    
    df_.columns = ["Start", "End"]
    df_["Chromosome"] = chromosome
    df_ = df_[["Chromosome", "Start", "End"]]

    features = seq_toolkit.GenomicFeatures(df_)
    features.merge()
    if out_path:
        features.write_bed(out_path)

    if return_df:
        return features.merged_df