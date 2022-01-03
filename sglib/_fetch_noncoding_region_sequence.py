import seq_toolkit

def _fetch_noncoding_region_sequence(regions_df, ref_genome_path):
    
    """"""
    
    chromosomes = regions_df.Chromosome.unique()    
    RegionSequences = {}
    ChromosomeSequences = {}
    
    for chromosome in chromosomes:
        chrom_df = regions_df.loc[regions_df.Chromosome == chromosome]
        chrom_start = chrom_df.Start.astype(int).min()
        chrom_end =  chrom_df.End.astype(int).max()
        
        ChromosomeSequences[chromosome] = seq_toolkit.fetch_chromosome(ref_genome_path, chromosome)
        RegionSequences[chromosome] = ChromosomeSequences[chromosome][chrom_start:chrom_end]
    
    del ChromosomeSequences[chromosome]
    
    return RegionSequences