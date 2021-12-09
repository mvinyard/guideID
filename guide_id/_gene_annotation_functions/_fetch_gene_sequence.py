
# _fetch_gene_sequence.py

__module_name__ = "_fetch_gene_sequence.py"
__author__ = ", ".join(["Michael E. Vinyard"])
__email__ = ", ".join(["vinyard@g.harvard.edu"])


# package imports #
# --------------- #
import seq_toolkit


# local imports #
# ------------- #
from ._return_gene_body_bounds import _return_gene_body_bounds
from ._fetch_chromosome_name import _fetch_chromosome_name


def _fetch_gene_sequence(ref_genome_path, gene_df):

    """
    Given a path to a reference genome fasta file and a gene_df, returns the sequence corresponding to that gene. 
    
    Parameters:
    -----------
    ref_genome_path
    
    gene_df
        obtained as a subset from gtfparse library and narrowed down to a single gene. 
    
    Returns:
    --------
    gene_seq
        gene sequence 
        type: str
    
    Notes:
    ------
    """

    gene_start, gene_end = _annotate_gene_body_start_stop(gene_df)
    chromosome = _get_chromosome_name(gene_df)
    chromosome_seq = seq_toolkit.fetch_chromosome(ref_genome_path, chromosome)
    gene_seq = chromosome_seq[gene_start:gene_end]

    return gene_seq