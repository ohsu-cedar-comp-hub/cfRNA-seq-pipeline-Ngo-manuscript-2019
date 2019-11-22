import pandas as pd


"""Accepts a HTSeq output directory and compiles {sample}_htseq_gene_count.txt into a joined tab sep file
Args:
    snakemake.input (list): list of globbed wildcards HTSeq 
    snakemake.output[0] (str): data/{params.project_id}_counts.txt
Returns:
    Compiled STAR gene counts table as tab delimited file.
"""

tables = [pd.read_csv(fh, sep='\t', index_col=0, names=[fh.split('/')[-1].split('_')[0]]) for fh in snakemake.input]
joined_table = pd.concat(tables, axis=1)
joined_sorted = joined_table.reindex(sorted(joined_table.columns), axis = 1)
joined_sorted.to_csv(snakemake.output[0], sep='\t')
