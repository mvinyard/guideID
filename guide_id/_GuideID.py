
# _GuideID.py

__module_name__ = "_GuideID.py"
__author__ = ", ".join(["Michael E. Vinyard"])
__email__ = ", ".join(["vinyard@g.harvard.edu",])


# local imports #
# ------------- #
from ._supporting_functions._annotate_protospacer_position import _annotate_protospacer_position
from ._supporting_functions._return_guides_in_regions import _return_guides_in_regions
from ._supporting_functions._FetchGene import _FetchGene


def _annotate_gene_body_start_stop(df):

    """Get start/stop locations"""
    gene_body = df.loc[df["feature"] == "gene"]
    return gene_body.start.values[0], gene_body.end.values[0]

class _GuideID:
    
    def __init__(self, fetch_gene=False, gene_feature="exon", reference_directory=False, save_feature_bed=True, df=False):
        
        """"""
        
        if fetch_gene:
            self.gene = _FetchGene(reference_directory)
            self.gene.fetch(fetch_gene, gene_feature, save_feature_bed)    
            self.sequence = self.gene.sequence
            self.guideID_df = self.gene.guideID_df
            self.gene_start, self.gene_end = _annotate_gene_body_start_stop(self.gene.gene_df)
        else:
            self.df = df
            self.global_start = df.sort_values("gene_feature.start")["gene_feature.start"].min()
            
        
    def scan(self,
             sequence=False,
             region_column="gene_feature",
             region_specification="exon",
             df=False,
             PAM="NGG",
             region_extension=0,
             return_guides=False,
            ):
        
        """
        Parameters:
        -----------
        df


        sequence


        region_column


        region_specification


        Returns:
        --------

        Notes:
        ------
        df.loc[df[region_column] == region_specification]
        """
        
        if df:
            self.df = df

        self.pam_df, self.guide_df, self.region_df = _return_guides_in_regions(
            self.sequence,
            self.guideID_df,
            region_column,
            region_specification,
            PAM,
            self.gene_start,
            region_extension,
        )
        
        self.guide_df = _annotate_protospacer_position(self.guide_df)
        
        if return_guides:
            return self.guide_df