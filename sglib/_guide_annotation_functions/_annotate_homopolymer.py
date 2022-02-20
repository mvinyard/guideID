
import numpy as np

def _annotate_homopolymer(df, protospacer_key="protospacer.sequence", n=4):

    """
    Annotate homopolymer sequences in a given protospacer.
    
    Parameters:
    -----------
    df
      sgRNA df
    
    protospacer_key
      key to the column containing the sgRNA protospacer sequences. 
    
    n
      Length of the homopolymer. 
      
    Returns:
    --------
    df_
      input df modified with homopolymer annotations.
    
    Notes:
    ------
    (1) Creates a new column in the DataFrame for each homopolymer sequence.
    (2) Currently, homopolymers are all passed at the same length. This could easily be changed in the future if needed. 
    """

    df_ = df.copy()

    for N in ["A", "C", "G", "T"]:
        df_['"poly{}".format(N)'] = df_[protospacer_key].str.contains(
            "".join(np.repeat(N, n))
        )

    return df_
