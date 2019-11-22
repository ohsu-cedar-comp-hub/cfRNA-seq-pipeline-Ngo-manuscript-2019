rule map_TE:
    input:
        fwd = "samples/trimmed/{sample}_R1_t.fq",
        rev = "samples/trimmed/{sample}_R2_t.fq"
    output:
        "samples/star_TE/{sample}/Aligned.out.bam",
        "samples/star_TE/{sample}/Log.final.out"
    params:
        index=config["star_index"],
        gtf=config["gtf_file"]
    run:
        STAR=config["star_tool"]

        shell("""
                {STAR} --runThreadN 12 --genomeDir {params.index} --sjdbGTFfile {params.gtf} \
                       --sjdbOverhang 100 --readFilesIn {input.fwd} {input.rev} \
                       --outFileNamePrefix samples/star_TE/{wildcards.sample}/ \
                       --outSAMtype BAM Unsorted --winAnchorMultimapNmax 200 --outFilterMultimapNmax 100""")

rule TEtranscripts:
    input:
        control = control_paths,
        treat = get_TE
    output:
        "results/TEtranscripts/{condition}.cntTable",
        "results/TEtranscripts/{condition}_sigdiff_gene_TE.txt"
    conda:
        "../envs/TE.yaml"
    params:
        gtf=config["gtf_file"],
        TE_gtf=config["TE_gtf"]
    shell:
        """TEtranscripts --format BAM --stranded reverse -t {input.treat} -c {input.control} \
                         --minread 1 -i 10 --padj 0.05 --GTF {params.gtf} --TE {params.TE_gtf} \
                         --mode multi --project results/TEtranscripts/{wildcards.condition}"""
