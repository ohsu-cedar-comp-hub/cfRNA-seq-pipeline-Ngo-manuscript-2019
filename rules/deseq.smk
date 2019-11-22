contrast = get_contrast

rule deseq2_init:
    input:
        counts = "data/{project_id}_counts.txt".format(project_id=config["project_id"])
    output:
        rds="results/diffexp/pairwise/{contrast}_all.rds",
        rld_out = "results/diffexp/pairwise/{contrast}_rlog_dds.rds"
    params:
        samples=config["omic_meta_data"],
        sample_id = config["sample_id"],
        linear_model = config["linear_model"],
        contrast = get_contrast
    conda:
        "../envs/deseq2.yaml"
    threads: get_deseq2_threads()
    script:
        "../scripts/deseq2-init.R"


rule deseq2_pairwise:
    input:
        rds="results/diffexp/pairwise/{contrast}_all.rds",
        rld = "results/diffexp/pairwise/{contrast}_rlog_dds.rds"
    output:
        table="results/diffexp/pairwise/{contrast}.diffexp.tsv",
        geneID_table="results/diffexp/pairwise/{contrast}.diffexp.geneID.tsv",
        ma_plot="results/diffexp/pairwise/{contrast}.ma_plot.pdf",
        p_hist="results/diffexp/pairwise/{contrast}.phist_plot.pdf",
        heatmap_plot = "results/diffexp/pairwise/{contrast}.heatmap_plot.pdf",
        panel_ma = "results/diffexp/pairwise/{contrast}.panel_ma.pdf",
        var_heat = "results/diffexp/pairwise/{contrast}.variance_heatmap.pdf",
        pca_plot = "results/diffexp/pairwise/{contrast}.pca_plot.pdf"
    params:
        contrast=get_contrast,
        linear_model = config["linear_model"],
        pca_labels = config["pca"]["labels"],
        sample_id = config["linear_model"]
    conda:
        "../envs/deseq2.yaml"
    threads: get_deseq2_threads()
    script:
        "../scripts/deseq2_pairwise.R"


rule deseq2_group:
    input:
        counts = "data/{project_id}_counts.txt".format(project_id = project_id)
    output:
        pca="results/diffexp/group/LRT_pca.pdf",
        sd_mean_plot="results/diffexp/group/LRT_sd_mean_plot.pdf",
        distance_plot = "results/diffexp/group/LRT_distance_plot.pdf",
        heatmap_plot = "results/diffexp/group/LRT_heatmap_plot.pdf",
        rds="results/diffexp/group/LRT_all.rds",
        rld_out = "results/diffexp/group/LRT_rlog_dds.rds"
    params:
        pca_labels=config["pca"]["labels"],
        samples=config["omic_meta_data"],
        sample_id = config["sample_id"],
        linear_model = config["linear_model"],
        LRT = config["diffexp"]["LRT"],
        colors = config['colors']['rcolorbrewer'],
        discrete = config['colors']['discrete']
    conda:
        "../envs/deseq2.yaml"
    script:
        "../scripts/deseq2_group.R"

rule deseq2_QC:
    input:
        rld="results/diffexp/group/LRT_rlog_dds.rds",
        rds="results/diffexp/group/LRT_all.rds"
    output:
        mds_plot="results/diffexp/group/MDS_plot.pdf",
        mds_table="results/diffexp/group/MDS_table.txt",
        heatmap_plot="results/diffexp/group/Heatmap_all_genes.pdf",
        sd_plot="results/diffexp/group/stdev_plot.pdf",
        rlogCounts_plot="results/diffexp/group/rlog_counts_violinPlot.pdf",
        rlogCounts_fac_plot="results/diffexp/group/rlog_counts_faceted_violinPlot.pdf",
        counts_plot="results/diffexp/group/counts_violinPlot.pdf",
        counts_fac_plot="results/diffexp/group/counts_faceted_violinPlot.pdf"
    params:
        sample_id = config["sample_id"],
        linear_model = config["linear_model"],
        colors = config['colors']['rcolorbrewer'],
        discrete = config['colors']['discrete']
    conda:
        "../envs/deseq2.yaml"
    script:
        "../scripts/QC.R"

rule deseq2_qplot:
    input:
        stats_table="results/diffexp/pairwise/{contrast}.diffexp.tsv",
    output:
        qplot="results/diffexp/pairwise/{contrast}.qplot.pdf",
        qhist="results/diffexp/pairwise/{contrast}.qhist.pdf",
        table = "results/diffexp/pairwise/{contrast}.qvalue_diffexp.tsv"
    params:
        contrast=get_contrast,
    conda:
        "../envs/omic_qc_wf.yaml"
    script:
        "../scripts/qplot.R"


rule deseq2_density:
    input:
        rld = "results/diffexp/group/LRT_rlog_dds.rds"
    output:
        density="results/diffexp/group/LRT_density_plot.pdf",
    params:
        linear_model = config["linear_model"],
        project_id = config["project_id"],
        colors = config['colors']['rcolorbrewer'],
        discrete = config['colors']['discrete']
    conda:
        "../envs/deseq2.yaml"
    script:
        "../scripts/density_plot.R"


rule GO:
    input:
        degFile="results/diffexp/pairwise/{contrast}.diffexp.tsv"
    output:
        "results/diffexp/pairwise/GOterms/{{contrast}}.diffexp.downFC.{FC}.adjp.{adjp}_BP_GO.txt".format(FC = config["FC"],adjp=config["adjp"]),
        "results/diffexp/pairwise/GOterms/{{contrast}}.diffexp.upFC.{FC}.adjp.{adjp}_BP_GO.txt".format(FC = config["FC"],adjp=config["adjp"])
    params:
        contrast = get_contrast,
        assembly = config["assembly"],
        printTree = config["printTree"],
        FC = config["FC"],
        adjp = config["adjp"]
    conda:
        "../envs/runGO.yaml"
    script:
        "../scripts/runGOforDESeq2.R"


rule volcano:
    input:
        degFile="results/diffexp/pairwise/{contrast}.diffexp.tsv"
    output:
        volcano_plot="results/diffexp/pairwise/{{contrast}}.diffexp.{adjp}.VolcanoPlot.pdf".format(adjp=config["adjp"])
    params:
        contrast = get_contrast,
        FC = config["FC"],
        adjp = config["adjp"]
    conda:
        "../envs/deseq2.yaml"
    script:
        "../scripts/RNAseq_makeVolcano.R"


rule permutation:
    input:
        counts = "data/{project_id}_counts.txt".format(project_id=config["project_id"])
    output:
        numGenes = "results/diffexp/pairwise/permutationTest/{contrast}.number.diff.genes.csv",
        permList = "results/diffexp/pairwise/permutationTest/{contrast}.permutation.list.csv",
        histogram = "results/diffexp/pairwise/permutationTest/Histogram.{contrast}.Permutation.Test.pdf"
    params:
        contrast = get_contrast,
        samples = config["omic_meta_data"],
        sample_id = config["sample_id"],
        linear_model = config["linear_model"]
    conda:
        "../envs/deseq2.yaml"
    script:
        "../scripts/permutation_test.R"


rule run_glimma:
    input:
        rds="results/diffexp/pairwise/{contrast}_all.rds"
    output:
        ma_plot = "results/diffexp/glimma-plots/{contrast}.ma_plot.html",
        volcano_plot = "results/diffexp/glimma-plots/{contrast}.volcano_plot.html",
    params:
        contrast = get_contrast,
        condition = config["linear_model"],
    conda:
        "../envs/glimma_env.yaml"
    script:
        "../scripts/run_glimma.R"


rule run_glimma_mds:
    input:
        rds="results/diffexp/group/LRT_all.rds"
    output:
        mds_plot = "results/diffexp/glimma-plots/{project_id}.mds_plot.html".format(project_id=project_id),
    params:
        project_id = config["project_id"],
    conda:
        "../envs/glimma_env.yaml"
    script:
        "../scripts/run_glimma_mds.R"
