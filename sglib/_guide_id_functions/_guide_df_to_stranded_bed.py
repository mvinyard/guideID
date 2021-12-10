def _write_single_strand_sgRNA_bed_file(df, gene, feature, strand, out_prefix):
    
    df_ = df.loc[df["pam.strand"] == strand]
    
    if strand == "+":
        strand_label = "forward"
    elif strand == "-":
        strand_label = "reverse"
    
    df_[["Chromosome", "protospacer.start", "protospacer.end"]].to_csv(
        "{}{}.{}.{}.bed".format(out_prefix, gene, feature, strand_label), header=False, index=False, sep="\t"
    )
    
def _guide_df_to_stranded_bed(df, gene, feature, out_prefix):
    
    """
    Save to two .bed files (top and bottom). 
    
    df
        guide_df
        
    path
        path prefix to output
    
    Notes:
    ------
    (1) It's useful to have separate strand bed files for visualization in IGV.
    """
    
    _write_single_strand_sgRNA_bed_file(df, gene, feature, "+", out_prefix)
    _write_single_strand_sgRNA_bed_file(df, gene, feature, "-", out_prefix)