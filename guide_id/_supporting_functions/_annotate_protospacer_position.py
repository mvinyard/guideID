
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

    """"""

    df_ = df.copy()
    df_ = _annotate_strand_multiplier(df, pam_strand="pam.strand")

    df_["protospacer.start"] = df_["pam.start"] + (20 * df_["strand.multiplier"])
    df_["protospacer.end"] = df_["pam.start"]
    df_ = df_.drop("strand.multiplier", axis=1)

    return df_