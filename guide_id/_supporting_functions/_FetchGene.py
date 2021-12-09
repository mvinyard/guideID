import gtfparse
import licorice
import os
import pandas as pd
import seq_toolkit
import warnings

warnings.filterwarnings("ignore")


def _check_gene_name(gene, gene_name_column):

    """See how many genes are returned from the GTF file based on the passed gene name."""

    if gene_name_column.nunique() > 1:

        gene_ = licorice.font_format(gene, ["BOLD", "CYAN"])

        print(
            "More than one gene name detected including the passed name: {}.\n".format(
                gene_
            )
        )
        for name in gene_name_column.unique():
            print("\t{}".format(gene_.join(name.split(gene))))
    else:
        return gene_name_column.unique()[0]


def _fetch_gene_annotation(gene, gtf):

    """Fetch gene annotation (feature boundaries) and the corresponding sequences."""

    gene_df = gtf.loc[gtf["gene_name"].str.contains(gene)]
    gene_id = _check_gene_name(gene, gene_df["gene_name"])

    return gene_df, gene_id


def _annotate_gene_body_start_stop(df):

    """Get start/stop locations"""
    gene_body = df.loc[df["feature"] == "gene"]
    return gene_body.start.values[0], gene_body.end.values[0]


def _merge_reduce_gene_region_bed(gene_df, feature_df, out_path=False, return_df=False):

    """"""

    chromosome = _get_chromosome_name(gene_df)
    df_ = feature_df.copy()    
    df_.columns = ["Start", "End"]
    df_["Chromosome"] = chromosome
    df_ = df_[["Chromosome", "Start", "End"]]

    features = seq_toolkit.GenomicFeatures(df_)
    features.merge()
    if out_path:
        features.write_bed(out_path)

    if return_df:
        return features.merged_df


def _isolate_feature(df, feature):

    """"""

    features = df.loc[df["feature"] == feature][["start", "end"]].values
    return pd.DataFrame(data=features, columns=["gene_feature.start", "gene_feature.end"])
    
def _get_chromosome_name(gene_df):

    """"""
    if gene_df.seqname.nunique() > 1:
        print("More than one chromosome identified:\n")
        for chromosome in gene_df.seqname.unique():
            print("\t{}".format(chromosome))
    else:
        return gene_df.seqname.unique()[0]


def _fetch_gene_sequence(ref_genome_path, gene_df):

    """"""

    gene_start, gene_end = _annotate_gene_body_start_stop(gene_df)
    chromosome = _get_chromosome_name(gene_df)
    chromosome_seq = seq_toolkit.fetch_chromosome(ref_genome_path, chromosome)
    gene_seq = chromosome_seq[gene_start:gene_end]

    return gene_seq

def _format_feature_df_for_guideID_search(feature_df, feature):

    df_ = feature_df.copy()
    start = df_.filter(regex="tart")
    end = df_.filter(regex="nd")
    df = pd.concat([start, end], axis=1)
    df.columns = ["gene_feature.start", "gene_feature.end"]

    df["gene_feature"] = feature
    df["length"] = df["gene_feature.end"] - df["gene_feature.start"]

    return df


class _FetchGene:
    
    """"""
    
    def __init__(self, reference_directory):

        """Assumes a 10x reference directory."""

        self.ref_genome_path = os.path.join(reference_directory, "fasta/genome.fa")
        self.gtf_path = os.path.join(reference_directory, "genes/genes.gtf")

        print("Loading GTF annotation file from {}...\n".format(self.gtf_path))
        self.gtf = gtfparse.read_gtf(self.gtf_path)

    def fetch(
        self,
        gene,
        feature,
        save_feature_bed=True,
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
        self.feature_df = _isolate_feature(self.gene_df, feature)
        print("Fetching gene sequence...\n")
        self.sequence = _fetch_gene_sequence(self.ref_genome_path, self.gene_df)
        if save_feature_bed:
            feature_bed_out_path = "{}.{}.bed".format(gene, feature)
            self.merged_feature_df = _merge_reduce_gene_region_bed(
                self.gene_df, self.feature_df, feature_bed_out_path, return_df=True
            )

        self.guideID_df = _format_feature_df_for_guideID_search(
            self.merged_feature_df, feature
        )