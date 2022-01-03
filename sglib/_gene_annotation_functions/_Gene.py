
# _Gene.py

__module_name__ = "_Gene.py"
__author__ = ", ".join(["Michael E. Vinyard"])
__email__ = ", ".join(["vinyard@g.harvard.edu"])


# import packages #
# --------------- #
import gtfparse
import os
import pandas as pd
import warnings


# local imports #
# ------------- #
from ._format_feature_df import _format_feature_df
from ._fetch_gene_annotation import _fetch_gene_annotation
from ._annotate_gene_feature import _annotate_gene_feature
from ._fetch_gene_sequence import _fetch_gene_sequence
from ._merge_reduce_gene_regions import _merge_reduce_gene_regions


# parameters #
# ---------- #
warnings.filterwarnings("ignore")


class _Gene:

    """"""

    def __init__(self, reference_gtf, reference_genome_path):

        """Assumes a 10x reference directory."""
        
        self.gtf = reference_gtf
        self.ref_genome_path = reference_genome_path
            
    def fetch(
        self,
        gene,
        feature,
        reverse_strand=False,
    ):

        """
        Get a DataFrame for a gene feature set.
        
        Parameters:
        -----------
        Gene
        
        Feature
        
        save_feature_bed
        
        """

        print("Fetching gene annotations...\n")
        self.gene_df, self.gene_id = _fetch_gene_annotation(gene, self.gtf)
        self.feature_df = _annotate_gene_feature(self.gene_df, feature)
        print("Fetching gene sequence...\n")
        self.sequence = _fetch_gene_sequence(self.ref_genome_path, self.gene_df)
        self.merged_feature_df = _merge_reduce_gene_regions(self.feature_df, gene, feature)
        self.formatted_feature_df = _format_feature_df(self.merged_feature_df, feature, reverse_strand)
        