import pandas as pd


"""Function accepts a STAR output directory and compiles all sample information from Log.final.out
Args:
    snakemake.input (list): list of globbed wildcards STAR Log.final.out
    project_title (str): Project title for compiled STAR mapping statistics

Returns:
    Compiled STAR log.final.out as tab delimited file.
"""

tables = [pd.read_csv(fh, sep = '\t', index_col = 0, names = [fh.split('/')[-2]]) for fh in snakemake.input]
joined_table = pd.concat(tables, axis=1)
joined_table_sorted = joined_table.reindex(sorted(joined_table.columns), axis = 1)
joined_table_sorted.to_csv(snakemake.output[0], sep='\t')
