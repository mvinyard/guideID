
def _annotate_GC_content(
    df, protospacer_key="protospacer.sequence", key_added="GC_content"
):

    """
    Parameters:
    -----------
    df
      sgRNA df

    protospacer_key
      key to the column containing the sgRNA protospacer sequences.

    key_added
      key indicating GC content.

    Returns:
    --------
    df_
      input df modified with GC content annotation.
    """

    df_ = df.copy()

    n_g = df_[protospacer_key].str.count("G")
    n_c = df_[protospacer_key].str.count("C")

    df_[key_added] = (n_g + n_c) * 100 / 20

    return df_
