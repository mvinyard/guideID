
from ._annotate_GC_content import _annotate_GC_content
from ._annotate_homopolymer import _annotate_homopolymer

def _annotate_biochemistry(
    df,
    protospacer_key="protospacer.sequence",
    GC_content_key_added="GC_content",
    homopolymer_n=4,
):

    """Adds information to the guide df about GC content and polynucleotide content"""

    df_ = df.copy()

    df_ = _annotate_homopolymer(df_, protospacer_key, key_added=GC_content_key_added)
    df_ = _annotate_GC_content(df_, protospacer_key, n=homopolymer_n)

    return df_
