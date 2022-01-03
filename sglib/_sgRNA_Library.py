
# _sgRNA_Library.py

__module_name__ = "_sgRNA_Library.py"
__author__ = ", ".join(["Michael E. Vinyard"])
__email__ = ", ".join(["vinyard@g.harvard.edu",])


# local imports #
# ------------- #
from ._gene_annotation_functions._Gene import _Gene
from ._gene_annotation_functions._return_gene_body_bounds import _return_gene_body_bounds
from ._gene_annotation_functions._print_exon_statistics import _print_exon_statistics

from ._guide_id_functions._return_guides_in_regions import _return_guides_in_regions
from ._guide_id_functions._annotate_protospacer_position import _annotate_protospacer_position
from ._guide_id_functions._guide_df_to_stranded_bed import _guide_df_to_stranded_bed
from ._guide_id_functions._plot_guide_distribution import _plot_guide_distribution

class _sgRNA_Library:
    
    def __init__(self, reference_directory):
        
        """"""
        
        self.reference_directory = reference_directory
            
    def from_gene(self, gene, gene_feature="exon", reverse_strand=False):
        
        """Build a set of regions based on an input gene name"""
        
        self.gene_name = gene
        self.gene_feature = gene_feature
        self.Gene = _Gene(self.reference_directory)
        self.Gene.fetch(gene, gene_feature, reverse_strand)
        self.sequence = self.Gene.sequence
        self.formatted_feature_df = self.Gene.formatted_feature_df
        self.gene_start, self.gene_end = _return_gene_body_bounds(self.Gene.gene_df)
        self.global_start = self.gene_start
        self.region_specification = gene_feature
        
        return _print_exon_statistics(self.formatted_feature_df)
        
    def from_regions(self, df, coordinate_keys=["Chromosome", "Start", "End"]):
        
        """
        Use a pre-built set of regions as input. 
        
        Parameters:
        -----------
        df
            type: pandas.DataFrame
            
        Returns:
        --------
        
        Notes:
        ------
        """
        
        self.formatted_df = df
        self.global_start = df.sort_values(coordinate_keys[1])[coordinate_keys[1]].min()
    
    def PAM_scan(self, PAM="NGG", extend_region=0, region_column="gene_feature", out_prefix="", return_guides=False,):
        
        """
        Parameters:
        -----------
        PAM
            CRISPR-Cas9 Protospacer Adjacent Motif sequence (though any motif would work!)
            default: "NGG"
            type: str
           
        extend_region
            number of nucleotides by which regions of interest should be extended. 
            type: int
            default: 0

        region_column
            prefix to the boundary indication columns. followed by ".start" and ".end"
            type: str
            default: "gene_feature"
    
        return_guides
            Boolean indicator if the resulting guide_df should be returned (in-memory subclass by default)
            default: False
            type: bool


        Returns:
        --------
        guide_df (or None by default)
        
        Notes:
        ------
        df.loc[df[region_column] == region_specification]
        """

        self.target_region_df = _return_guides_in_regions(
            self.sequence,
            self.formatted_feature_df,
            region_column,
            self.region_specification,
            PAM,
            self.global_start,
            extend_region,
        )         
        
        self.guide_df = _annotate_protospacer_position(self.target_region_df)
        _guide_df_to_stranded_bed(self.guide_df, self.gene_name, self.gene_feature, out_prefix)
        
        _plot_guide_distribution(self.guide_df)
        self.guide_df.to_excel("./{}.{}.{}.sgRNAs.xlsx".format(self.gene_name, self.gene_feature, PAM), index=False)
        
        if return_guides:
            return self.guide_df