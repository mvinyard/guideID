
from gtfparse import read_gtf
import os
import pandas as pd

def _parse_gtf(reference_directory):
    
    """
    Given a path to a reference directory with structure:
    
        /path/to/reference_directory/fasta/genome.fa
        /path/to/reference_directory/genes/genes.gtf
    
    Parameters:
    -----------
    reference_directory
    
    Returns:
    --------
    gtf
    
    Creates: /path/to/reference_directory/genes/gtf.tsv
    
    Notes:
    ------
    
    """

    ref_genome_path = os.path.join(reference_directory, "fasta/genome.fa")
    gtf_path = os.path.join(reference_directory, "genes/genes.gtf")
    gtf_tsv = os.path.join(reference_directory, "genes/gtf.tsv")

    if os.path.exists(gtf_tsv):
        print("Loading GTF annotation file from {}...\n".format(gtf_tsv))
        gtf = pd.read_csv(gtf_tsv, sep="\t")
    else:
        print("Loading GTF annotation file from {}...\n".format(gtf_path))
        gtf = read_gtf(gtf_path)
        gtf[["seqname",
                 "feature",
                 "gene_type",
                 "gene_name",
                 "start",
                 "end",
                 "strand",
                 "exon_number",]].to_csv(gtf_tsv, sep="\t", index=False)
        
    return gtf