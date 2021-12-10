
__module_name__ = "_fetch_gene_annotation.py"
__author__ = ", ".join(["Michael E. Vinyard"])
__email__ = ", ".join(["vinyard@g.harvard.edu"])


# import packages #
# --------------- #
import licorice


def _check_gene_name(gene, gene_name_column):

    """
    See how many genes are returned from the GTF file based on the passed gene name.
        
    Parameters:
    -----------
    if only one gene name is found (success): gene_id
    
    if >1 gene name is found: None
        Instead prints the alternative gene names found, highlighting where the passed gene name overlaps non-uniquely.
    """

    if gene_name_column.nunique() > 1:

        gene_ = licorice.font_format(gene, ["BOLD", "CYAN"])

        print(
            "More than one gene name detected including the passed name: {}.\n".format(
                gene_
            )
        )
        for name in gene_name_column.unique():
            print("\t{}".format(gene_.join(name.split(gene))))
    else:
        return gene_name_column.unique()[0]
    
    
def _fetch_gene_annotation(gene, gtf):

    """
    Fetch gene annotation (feature boundaries) and the corresponding sequences.
    
    Parameters:
    -----------
    gene
        gene name that should be found in the "gene_name" column of the GTF DataFrame.
        type: str
        
    gtf
        GTF annotation DataFrame loaded by the gtfparse library. 
        pandas.DataFrame
    
    Returns:
    --------
    gene_df
        subset of the input gtf DataFrame corresponding to rows that match the input gene
        type: pandas.DataFrame
    
    gene_id
        name of the gene. ideally mathces the passed "gene" argument. 
        type: str
    """

    gene_df = gtf.loc[gtf["gene_name"].str.contains(gene)]
    gene_id = _check_gene_name(gene, gene_df["gene_name"])

    return gene_df, gene_id