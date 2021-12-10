
import seq_toolkit

def _merge_reduce_gene_regions(feature_df, gene, feature):

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
    """
    feature_bed_out_path = "{}.{}.bed".format(gene, feature)
    
    df_ = feature_df.copy()    
    df_.columns = ["Start", "End", "Chromosome"]

    features = seq_toolkit.GenomicFeatures(df_[["Chromosome", "Start", "End"]])
    features.merge()
    features.write_bed(feature_bed_out_path)
    
    return features.merged_df