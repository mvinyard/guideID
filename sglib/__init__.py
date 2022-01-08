# __init__.py

from ._sgRNA_Library import _sgRNA_Library as Library

# functions
from ._gene_annotation_functions._merge_reduce_gene_regions import _merge_reduce_gene_regions as merge_reduce

from ._guide_annotation_functions._annotate_biochemistry import _annotate_biochemistry as annotate_biochemistry, _annotate_GC_content as annotate_GC, _annotate_homopolymer as annotate_homopolymer
