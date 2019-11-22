rule trimming:
    input:
        fwd = "samples/raw/{sample}_R1.fq",
        rev = "samples/raw/{sample}_R2.fq"
    output:
        fwd = "samples/trimmed/{sample}_R1_t.fq",
        rev = "samples/trimmed/{sample}_R2_t.fq",
        single = "samples/trimmed/{sample}_R1_singletons.fq"
    run:
        sickle = config["sickle_tool"]

        shell("{sickle} pe -f {input.fwd} -r {input.rev}  -l 40 -q 20 -t sanger  -o {output.fwd} -p {output.rev} -s {output.single} &> {input.fwd}.log")

rule fastqc:
    input:
        fwd = "samples/trimmed/{sample}_R1_t.fq",
        rev = "samples/trimmed/{sample}_R2_t.fq"
    output:
        fwd = "samples/fastqc/{sample}/{sample}_R1_t_fastqc.zip",
        rev = "samples/fastqc/{sample}/{sample}_R2_t_fastqc.zip"
    conda:
        "../envs/fastqc.yaml"
    shell:
        """fastqc --outdir samples/fastqc/{wildcards.sample} --extract  -f fastq {input.fwd} {input.rev}"""

rule fastqscreen:
    input:
        fwd = "samples/trimmed/{sample}_R1_t.fq",
        rev = "samples/trimmed/{sample}_R2_t.fq"
    output:
        "samples/fastqscreen/{sample}/{sample}_R1_t_screen.html",
        "samples/fastqscreen/{sample}/{sample}_R1_t_screen.png",
        "samples/fastqscreen/{sample}/{sample}_R1_t_screen.txt",
        "samples/fastqscreen/{sample}/{sample}_R2_t_screen.html",
        "samples/fastqscreen/{sample}/{sample}_R2_t_screen.png",
        "samples/fastqscreen/{sample}/{sample}_R2_t_screen.txt"
    params:
        conf = config["conf"]
    conda:
        "../envs/fastqscreen.yaml"
    shell:
        """fastq_screen --aligner bowtie2 --conf {params.conf} --outdir samples/fastqscreen/{wildcards.sample} {input.fwd} {input.rev}"""

rule star:
    input:
        fwd = "samples/trimmed/{sample}_R1_t.fq",
        rev = "samples/trimmed/{sample}_R2_t.fq"
    output:
        "samples/star/{sample}_bam/Aligned.sortedByCoord.out.bam",
        "samples/star/{sample}_bam/ReadsPerGene.out.tab",
        "samples/star/{sample}_bam/Log.final.out"
    params:
        gtf=config["gtf_file"]
    run:
         STAR=config["star_tool"],
         pathToGenomeIndex = config["star_index"]

         shell("""
                {STAR} --runThreadN {threads} --runMode alignReads --genomeDir {pathToGenomeIndex} \
                --readFilesIn {input.fwd} {input.rev} \
                --outFileNamePrefix samples/star/{wildcards.sample}_bam/ \
                --sjdbGTFfile {params.gtf} --quantMode GeneCounts \
                --sjdbGTFtagExonParentGene gene_name \
                --outSAMtype BAM SortedByCoordinate \
                #--readFilesCommand zcat \
                --twopassMode Basic
                """)

rule index:
    input:
        "samples/star/{sample}_bam/Aligned.sortedByCoord.out.bam"
    output:
        "samples/star/{sample}_bam/Aligned.sortedByCoord.out.bam.bai"
    conda:
        "../envs/omic_qc_wf.yaml"
    shell:
        """samtools index {input} {output}"""

rule star_statistics:
    input:
        expand("samples/star/{sample}_bam/Log.final.out",sample=SAMPLES)
    output:
        "results/tables/{project_id}_STAR_mapping_statistics.txt".format(project_id = config["project_id"])
    script:
        "../scripts/compile_star_log.py"


rule picard:
    input:
        "samples/star/{sample}_bam/Aligned.sortedByCoord.out.bam"
    output:
        temp("samples/genecounts_rmdp/{sample}_bam/{sample}.rmd.bam")
    params:
        name="rmd_{sample}",
        mem="5300"
    run:
      picard=config["picard_tool"]

      shell("java -Xmx3g -jar {picard} \
      INPUT={input} \
      OUTPUT={output} \
      METRICS_FILE=samples/genecounts_rmdp/{wildcards.sample}_bam/{wildcards.sample}.rmd.metrics.text \
      REMOVE_DUPLICATES=true")


rule sort:
    input:
      "samples/genecounts_rmdp/{sample}_bam/{sample}.rmd.bam"
    output:
      "samples/genecounts_rmdp/{sample}_bam/{sample}_sort.rmd.bam"
    params:
      name = "sort_{sample}",
      mem = "6400"
    conda:
      "../envs/omic_qc_wf.yaml"
    shell:
      """samtools sort -O bam -n {input} -o {output}"""


rule samtools_stats:
    input:
        "samples/genecounts_rmdp/{sample}_bam/{sample}_sort.rmd.bam"
    output:
        "samples/samtools_stats/{sample}.txt"
    log:
        "logs/samtools_stats/{sample}_samtools_stats.log"
    conda:
        "../envs/omic_qc_wf.yaml"
    wrapper:
        "0.17.0/bio/samtools/stats"


rule genecount:
    input:
        "samples/genecounts_rmdp/{sample}_bam/{sample}_sort.rmd.bam"
    output:
        "samples/htseq_count/{sample}_htseq_gene_count.txt",
    log:
        "logs/genecount/{sample}_genecount.log"
    params:
        name = "genecount_{sample}",
        gtf = config["gtf_file"]
    conda:
        "../envs/omic_qc_wf.yaml"
    threads: 1
    shell:
        """htseq-count \
                -f bam \
                -r name \
                -s reverse \
                -m intersection-strict \
                {input} \
                {params.gtf} > {output}"""


rule count_exons:
    input:
        "samples/genecounts_rmdp/{sample}_bam/{sample}_sort.rmd.bam"
    output:
        "samples/htseq_exon_count/{sample}_htseq_exon_count.txt"
    params:
        exon_gtf = config["exon_gtf"]
    conda:
        "../envs/omic_qc_wf.yaml"
    shell:
        """htseq-count \
               -f bam \
               -r name \
               -s reverse \
               -m intersection-strict \
               -i exon_id \
               --additional-attr=gene_name \
               {input} \
               {params.exon_gtf} > {output}"""


rule compile_counts:
    input:
        expand("samples/htseq_count/{sample}_htseq_gene_count.txt",sample=SAMPLES)
    output:
        "data/{project_id}_counts.txt".format(project_id=config["project_id"])
    script:
        "../scripts/compile_counts_table.py"


rule compile_counts_and_stats:
    input:
        expand("samples/htseq_count/{sample}_htseq_gene_count.txt",sample=SAMPLES)
    output:
        "data/{project_id}_counts_w_stats.txt".format(project_id=config["project_id"])
    script:
        "../scripts/compile_counts_table_w_stats.py"


rule compile_exon_counts:
    input:
        expand("samples/htseq_exon_count/{sample}_htseq_exon_count.txt", sample=SAMPLES)
    output:
        "data/{project_id}_exon_counts.txt".format(project_id = config["project_id"])
    script:
        "../scripts/compile_exon_counts.R"

