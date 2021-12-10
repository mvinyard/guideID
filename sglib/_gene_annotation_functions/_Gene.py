
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

    def __init__(self, reference_directory):

        """Assumes a 10x reference directory."""

        self.ref_genome_path = os.path.join(reference_directory, "fasta/genome.fa")
        self.gtf_path = os.path.join(reference_directory, "genes/genes.gtf")
        self.gtf_tsv = os.path.join(reference_directory, "genes/gtf.tsv")

        if os.path.exists(self.gtf_tsv):
            print("Loading GTF annotation file from {}...\n".format(self.gtf_tsv))
            self.gtf = pd.read_csv(self.gtf_tsv, sep="\t")
        else:
            print("Loading GTF annotation file from {}...\n".format(self.gtf_path))
            self.gtf = gtfparse.read_gtf(self.gtf_path)
            self.gtf[["seqname",
                     "feature",
                     "gene_type",
                     "gene_name",
                     "start",
                     "end",
                     "strand",
                     "exon_number",]].to_csv(self.gtf_tsv, sep="\t", index=False)
            
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
        