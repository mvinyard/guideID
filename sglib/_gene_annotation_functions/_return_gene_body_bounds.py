
__module_name__ = "_return_gene_body_bounds.py"
__author__ = ", ".join(["Michael E. Vinyard"])
__email__ = ", ".join(["vinyard@g.harvard.edu"])

def _return_gene_body_bounds(df):

    """
    Get start/end locations of a gene
    
    Parameters:
    -----------
    df
        df that is subset from reading a GTF file using the GTFparse library.
    
    Returns:
    --------
    start, end
    
    Notes:
    ------
    (1) This function could be extended to other features though would require an additioanl level
        of identifying the max/min of feature sets that have multiple rows in a GTF (as opposed to 
        the "gene" feature, which has only a single row / set of values. 
    """
    
    gene_body = df.loc[df["feature"] == "gene"]
    
    return gene_body.start.values[0], gene_body.end.values[0]