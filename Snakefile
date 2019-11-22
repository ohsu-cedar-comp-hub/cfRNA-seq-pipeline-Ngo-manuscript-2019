__author__ = "Joey Estabrook"
__email__ = "estabroj@ohsu.edu"
__license__ = "MIT"

"""Computation Hub omic data processing pipeline"""

import datetime
import sys
import os
import pandas as pd
import json

timestamp = ('{:%Y-%m-%d_%H:%M:%S}'.format(datetime.datetime.now()))

configfile:"omic_config.yaml"
project_id = config["project_id"]

SAMPLES, = glob_wildcards("samples/raw/{sample}_R1.fq")

md = pd.read_table(config["omic_meta_data"], index_col="PP_ID",dtype=str)
condition = config["linear_model"]
baseline = config["TE_baseline"]

# TO FILTER 
control = md[md[condition]==baseline]
treat = md.loc[~md.index.isin(control.index)] 

control_paths = ['samples/star_TE/{}/Aligned.out.bam'.format(x) for x in control.index]

treat['cond_paths'] = ['samples/star_TE/{}/Aligned.out.bam'.format(x) for x in treat.index]
cond_paths = {key:x['cond_paths'].values.tolist() for key,x in treat.groupby(condition)}
CONDITIONS = cond_paths.keys()

# Wildcard function to grab proper condition
def get_TE(wildcards):
    return cond_paths[wildcards.condition]

ext = ['r','R1.pdf','R2.pdf','xls']
fastq_ext = ['R1','R2']
fastqscreen_ext = ['html','png','txt']
insertion_and_clipping_prof_ext = ['r','R1.pdf','R2.pdf','xls']
inner_distance_ext = ['_freq.txt','_plot.pdf','_plot.r','.txt']
read_dist_ext = ['txt']
read_gc_ext = ['.xls','_plot.r','_plot.pdf']


with open('cluster.json') as json_file:
    json_dict = json.load(json_file)

rule_dirs = list(json_dict.keys())

for rule in rule_dirs:
    if not os.path.exists(os.path.join(os.getcwd(),'logs',rule)):
        log_out = os.path.join(os.getcwd(), 'logs', rule)
        os.makedirs(log_out)
        print(log_out)


result_dirs = ['diffexp','tables']
for rule in result_dirs:
    if not os.path.exists(os.path.join(os.getcwd(),'results',rule)):
        log_out = os.path.join(os.getcwd(), 'results', rule)
        os.makedirs(log_out)
        print(log_out)

def message(mes):
    sys.stderr.write("|--- " + mes + "\n")


def get_deseq2_threads(wildcards=None):
    few_coeffs = False if wildcards is None else len(get_contrast(wildcards)) < 10
    return 1 if len(config["omic_meta_data"]) < 100 or few_coeffs else 6


def get_contrast(wildcards):
    """Return each contrast provided in the configuration file"""
    return config["diffexp"]["contrasts"][wildcards.contrast]


for sample in SAMPLES:
    message("Sample " + sample + " will be processed")

for condition in CONDITIONS:
    message("Condition " + condition + " will be processed")
#print(control_paths)
#print(cond_paths)

rule all:
    input:
        expand("results/tables/{project_id}_STAR_mapping_statistics.txt", project_id = config['project_id']),
        expand("samples/fastqc/{sample}/{sample}_{fastq_ext}_t_fastqc.zip", sample = SAMPLES, fastq_ext = fastq_ext),
        expand("samples/fastqscreen/{sample}/{sample}_{fastq_ext}_t_screen.{fastqscreen_ext}", sample=SAMPLES, fastq_ext=fastq_ext, fastqscreen_ext=fastqscreen_ext),
        "data/{project_id}_counts_w_stats.txt".format(project_id=config['project_id']),
        "data/{project_id}_exon_counts.txt".format(project_id = config["project_id"]),
        expand("rseqc/insertion_profile/{sample}/{sample}.insertion_profile.{ext}",sample=SAMPLES, ext=insertion_and_clipping_prof_ext),
        expand("rseqc/inner_distance/{sample}/{sample}.inner_distance{ext}", sample = SAMPLES, ext = inner_distance_ext),
        expand("rseqc/clipping_profile/{sample}/{sample}.clipping_profile.{ext}", sample = SAMPLES, ext = insertion_and_clipping_prof_ext),
        expand("rseqc/read_distribution/{sample}/{sample}.read_distribution.{ext}", sample = SAMPLES, ext = read_dist_ext),
        expand("rseqc/read_GC/{sample}/{sample}.GC{ext}", sample = SAMPLES, ext = read_gc_ext),
        "results/tables/read_coverage.txt",
        expand("samples/star_TE/{sample}/Aligned.out.bam", sample = SAMPLES),
        expand("results/TEtranscripts/{condition}.cntTable", condition = CONDITIONS),
        expand("results/diffexp/pairwise/{contrast}.pca_plot.pdf", contrast = config["diffexp"]["contrasts"]),
        "results/diffexp/group/LRT_pca.pdf",
        "results/diffexp/group/MDS_table.txt",
        "results/diffexp/group/LRT_density_plot.pdf",
        expand(["results/diffexp/pairwise/{contrast}.qplot.pdf","results/diffexp/pairwise/{contrast}.qhist.pdf","results/diffexp/pairwise/{contrast}.qvalue_diffexp.tsv"],contrast=config["diffexp"]["contrasts"]),
        expand(["results/diffexp/pairwise/GOterms/{contrast}.diffexp.downFC.{FC}.adjp.{adjp}_BP_GO.txt", "results/diffexp/pairwise/GOterms/{contrast}.diffexp.upFC.{FC}.adjp.{adjp}_BP_GO.txt"], contrast = config["diffexp"]["contrasts"], FC=config['FC'], adjp=config['adjp']),
        expand("results/diffexp/pairwise/{contrast}.diffexp.{adjp}.VolcanoPlot.pdf", contrast = config["diffexp"]["contrasts"], adjp = config['adjp']),
        expand("results/diffexp/pairwise/permutationTest/Histogram.{contrast}.Permutation.Test.pdf", contrast = config["diffexp"]["contrasts"]),
        expand(["results/diffexp/glimma-plots/{contrast}.ma_plot.html", "results/diffexp/glimma-plots/{contrast}.volcano_plot.html"],contrast = config["diffexp"]["contrasts"]),
        "results/diffexp/glimma-plots/{project_id}.mds_plot.html".format(project_id=project_id)

include: "rules/align_rmdp.smk"
include: "rules/omic_qc.smk"
include: "rules/TE.smk"
include: "rules/deseq.smk"
