
rule transcript_abundance_plots:
    input:
        expand("results/tmp/{sample.sample}.ready.for.plots", sample=samples.reset_index().itertuples())
    output:
        txi="results/txi.RData",
        heatmap="results/Heatmap_Most_Var.png",
        pca="results/PCA.png"
    params:
        current_dir=get_cwd(),
        txgenes=resolve_single_filepath(*references_abs_path(ref="references"), config.get("tx2gene_hsa"))
    conda:
        "../envs/rplots.yaml"
    script:
        "scripts/Analysis_v1.R"
