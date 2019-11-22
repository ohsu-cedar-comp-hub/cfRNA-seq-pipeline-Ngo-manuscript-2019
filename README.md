cfRNA-seq Pipeline
======================

Pipeline to assess quality of Omics datasets, specifically tailored towards cell-free RNA.

This is a package of Python and R scripts that enable reading, processing and analysis of cfRNA Omics' datasets. 
This package implements the Snakemake management workflow system and is currently implemented to work with 
the cluster management and job scheduling system SLURM. This snakemake workflow utilizes conda installations to download and use packages for further analysis, so please ensure that you have installed miniconda prior to use.

Questions/issues
======================

Please add an issue to the cfRNA-seq repository. We would appreciate if your issue included sample code/files 
(as appropriate) so that we can reproduce your bug/issue. 


Contributing
======================

We welcome contributors! For your pull requests, please include the following:

* Sample code/file that reproducibly causes the bug/issue
* Documented code providing fix
* Unit tests evaluating added/modified methods. 

Use
======================

Locate raw files:
* After sequencing, your raw fastq files are placed in `/path/to/sequencing/files`.

```
$ cd /path/to/raw/data
$ ls -alh
```

Check md5sum.

```
$ md5sum –c md5sum.txt > md5sum_out.txt
```

Move your files into the archive to be stored.

```
$ mv /path/to/raw/data /path/to/archive
```

Check md5sum again to ensure your sequencing files are not corrupted.

```
$ md5sum –c md5sum.txt > md5sum_out.txt
```

Unzip all fastq files.

```
$ gunzip –d sample.fastq.gz
$ ctrl+z
$ bg
```

Clone the Omics-QC Pipeline into your working directory.

```
$ git clone https://github.com/ohsu-cedar-comp-hub/Omics-QC-pipeline.git
```

Create a `samples/raw` directory, a `logs` directory and a `data` directory (if they do not exist) in your `wdir()`.

```
$ mkdir logs
$ mkdir data
$ mkdir samples
$ cd samples
$ mkdir raw
```

Symbollically link the fastq files of your samples to the `wdir/samples/raw` directory using a bash script loop in your terminal.

```
$ ls -1 /path/to/data/LIB*R1*fastq | while read fastq; do
    R1=$( basename $fastq | cut -d _ -f 2 | awk '{print $1"_R1.fq"}' )
    R2=$( basename $fastq | cut -d _ -f 2 | awk '{print $1"_R2.fq"}' )
    echo $R1 : $R2
    ln -s $fastq ./$R1
    ln -s ${fastq%R1_001.fastq}R2_001.fastq ./$R2
done
```

Upload your metadata file to the `data` directory, with the correct formatting:
* Columns should read:
```StudyID  Column2   Column3   ...```
* Each row should be a sample, with subsequent desired information provided (RNA extraction date, etc.)
* All values in this file should be tab-separated

Edit the `omic_config.yaml` in your `wdir()`:
* Change the `project_id` to a unique project identifier
* Add appropriate contrasts based on your samples under the `[diffexp][contrasts]` section
* Add the path to your metadata file for the `omic_meta_data` and `samples` parameters
* Change `base_dir` to your current working directory
* Ensure you have the correct `assembly` specified
    * Current options for this are: hg19, hg38.89 (ensembl v89) and hg38.90 (ensembl v90)

Do a dry-run of snakemake to ensure proper execution before submitting it to the cluster (in your wdir).

```
$ snakemake -np --verbose
```

Once your files are symbolically linked, you can submit the job to exacloud via your terminal window.

```
$ sbatch submit_snakemake.sh
```

To see how the job is running, look at your queue.

```
$ squeue -u your_username
```

Detailed Workflow
=================================

Alignment
======================
1) Trimming
    * Trimming of paired-end reads was performed using the trimming tool `sickle`
    * The output is located in `samples/trimmed/`
2) Quality Analysis
    * Trimmed reads were subject to `fastqc` quality analysis
    * The output is located in `samples/fastqc/{sample}/{samples}_t_fastqc.zip`
3) Alignment
    * Trimmed reads were aligned to the hg38 genome assembly using `STAR`
        * We included a two pass mode flag in order to increase the number of aligned reads
        * Output is placed in `samples/star/{sample}_bam/`
            * Output directory includes: `Aligned.sortedByCoord.out.bam`, `ReadsPerGene.out.tab`, and `Log.final.out`
    * We extracted the statistics from the `STAR` run, and placed them in a table, summarizing the results across all samples from the `Log.final.out` output of STAR
        * Output is `results/tables/{project_id}_STAR_mapping_statistics.txt`
4) Summarizing output
    * `htseq` and `samtools` were used to extract the gene counts for each sample, and `picard` was used to remove duplicate reads
    * We summarize these results into 1 table, which includes the gene counts across all samples
    * The output is located in `data/{project_id}_counts.txt`

Quality Analysis / Quality Check
======================
1) RSEQC Quality check 
    * `RSEQC` was used to check the quality of the reads by using a collection of commands from the `RSEQC` package:
        * Insertion Profile
        * Inner Distance
        * Clipping Profile
        * Read distribution
        * Read GC
    * For more information on these, visit: http://dldcc-web.brc.bcm.edu/lilab/liguow/CGI/rseqc/_build/html/index.html#usage-information
    * Output directory: `rseqc/`
2) QA/QC scripts to analyze the data as a whole 
    * The purpose of this analysis is to identify potential batch effects and outliers in the data
    * The outputs to this are located in the `results` directory, and are distributed amongst 4 subdirectories, numbered `1 through 4`
        * `1`
            * A *boxplot* of the raw log2-transformed gene counts across all samples
            * A *boxplot* of the loess-transformed gene counts across all samples
            * A *scatter plot* comparing raw gene counts to loess-transformed gene counts
            * A *density plot* of raw log2-transformed gene counts across all samples 
            * A *density plot* of loess-transformed gene counts across all samples
            * A *scatter plot* of the standard deviation of raw log2-transformed gene counts across all samples
            * A *scatter plot* of the standard deviation of loess-transformed gene counts across all samples
        * `2`
            * A *heatmap* of all raw log2-transformed gene counts across samples
            * A *heatmap* of all loess-transformed gene counts across samples
                * These are generated to look for any batch effects in the data, due to date of extraction, or other factors
            * An *MDS Plot* for all samples, generated with the raw log2-transformed gene counts
            * An *MDS Plot* for all samples, generated with the loess-transformed gene counts
                * These are generated to look for outliers in the data
        * `3`
            * *p-value histograms* for each contrast specified in the `omic_config.yaml`
            * *q-value QC plot arrays* for each contrast specified in the `omic_config.yaml`
        * `4`
            * A *Heatmap* which looks at genes with a high FC and low q-value (very significant)
                * Takes genes with a FC>1.3, and ranks those by q-value. From this, a heatmap is generated for the top *50, 100 and 200* genes in this list
            * An *MDS Plot* which looks at the same subsets of genes as the Heatmap described above
            
Differential Expression Analysis (DESeq2)
======================
1) Initializing the DESeq2 object
    * Here, we run `DESeq2` on the genecounts table, which generates an RDS object and rlog
        * This includes the DE analysis across all samples
        * Output is located in the `results/diffexp/ directory`
    * From the dds object generated, we extract the normalized counts and generate a table with the results
        * Output is `results/tables/{project_id}_normed_counts.txt`
2) Generating plots
    * From the RDS object, we generate a collection of informative plots. These include:
        * *PCA Plot*
        * *Standard Deviation from the Mean Plot*
        * *Heatmap*
        * *Variance Heatmap*
        * *Distance Plot*
3) Differential Expression Analysis
    * We perform Differential Expression (DE) analysis for each contrast listed in the `omic_config.yaml`
    * Our output consists of DE gene count tables and a variety of plots
        * A table is generated for genes that are differentially expressed for each contrast
            * The output is placed in `results/diffexp/{contrast}.diffexp.tsv`
        * *MA Plots* are generated for each contrast
        * *p-histograms* are generated for each contrast
4) Differential Expression Plots
    * We use the output from DESeq2 to generate two types of plots:
        * Gene Ontology (GO) plots:
            * A `tree graph` describing the GO ID relationship for significantly up/downregulated genes in a given comparison
                * Output is located in `results/diffexp/GOterms`
            * A `bar graph` describing the enrichment and significance of GO IDs for up/downregulated genes in a given comparison
        * Volcano plots:
            * A `volcano plot` describing the distribution of up/downregulated genes in a given comparison
                * Output is located in `results/diffexp`
