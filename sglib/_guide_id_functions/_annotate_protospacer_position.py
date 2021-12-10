
__module_name__ = "_annotate_protospacer_position.py"
__author__ = ", ".join(["Michael E. Vinyard"])
__email__ = ", ".join(["vinyard@g.harvard.edu"])


# import packages #
# --------------- #
import pandas as pd


def _strand_multiplier_df(pam_strand="pam.strand"):

    """"""

    strand_multiplier = pd.DataFrame(
        data=[["+", -1], ["-", 1]],
        columns=[
            pam_strand,
            "strand.multiplier",
        ],
    )

    return strand_multiplier


def _annotate_strand_multiplier(df, pam_strand="pam.strand"):

    df_ = df.copy()
    strand_multiplier = _strand_multiplier_df(pam_strand)
    return df_.merge(strand_multiplier, on=pam_strand)


def _annotate_protospacer_position(df):

    """
    Given a pandas DataFrame as input with "pam.start" and "pam.strand" columns, protospacer 
    position (not sequence) is then annotated.
    
    Parameters:
    -----------
    df
        type: pandas.DataFrame
    
    Returns:
    --------
    df
        Modified version of the input pandas DataFrame with the added columns: ["protospacer.start", "protospacer.end"]
        type: pandas.DataFrame
        
    Notes:
    ------
    (1) The strand.multiplier column is added and dropped before the final return of the modified DataFrame. 
    
    """

    df_ = df.copy()
    df_ = _annotate_strand_multiplier(df, pam_strand="pam.strand")

    df_["protospacer.start"] = df_["pam.start"] + (20 * df_["strand.multiplier"])
    df_["protospacer.end"] = df_["pam.start"]
    df_ = df_.drop("strand.multiplier", axis=1)

    return df_