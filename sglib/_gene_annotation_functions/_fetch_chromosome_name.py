
def _fetch_chromosome_name(gene_df):

    """"""
    
    if gene_df.seqname.nunique() > 1:
        print("More than one chromosome identified:\n")
        for chromosome in gene_df.seqname.unique():
            print("\t{}".format(chromosome))
    else:
        return gene_df.seqname.unique()[0]
