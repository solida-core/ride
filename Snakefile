import pandas as pd
from snakemake.utils import validate, min_version
##### set minimum snakemake version #####
min_version("5.1.2")


##### load config and sample sheets #####
#configfile: "config.yaml"

## USER FILES ##
samples = pd.read_csv(config["samples"], index_col="sample", sep="\t")
## ---------- ##


include:
    "rules/functions.py"

rule all:
    input:
        # Kallisto+Sleuth workflow
#        "kallisto/DEGS/gene_table.txt",
        expand("kallisto/{sample.sample}/abundance.tsv", sample=samples.reset_index().itertuples()),
        # Star
        expand("star/{sample.sample}/{sample.sample}.Aligned.sortedByCoord.out.bam", sample=samples.reset_index().itertuples()),
        # Rseqc
        expand("rseqc/{sample.sample}/{sample.sample}.saturation.pdf", sample=samples.reset_index().itertuples()),
        "qc/multiqc.html",
        expand("star/{sample.sample}/count/{sample.sample}_featurecounts.cnt",sample=samples.reset_index().itertuples()),
        expand("star/{sample.sample}/count/{sample.sample}_HTSeqcounts.cnt",sample=samples.reset_index().itertuples()),
        "results/Heatmap_Most_Var.png"

include_prefix="rules"
include:
    include_prefix + "/rseqc.smk"
include:
    include_prefix + "/plots.smk"
if config.get("read_type")=="se":
    include:
        include_prefix + "/trimming_se.smk"
    include:
        include_prefix + "/qc_se.smk"
    include:
        include_prefix + "/reads_feature_count.smk"
    include:
        include_prefix + "/kallisto.smk"
    include:
        include_prefix + "/star2.smk"

else:
    include:
        include_prefix + "/pe/trimming_pe.smk"
    include:
        include_prefix + "/pe/qc_pe.smk"
    include:
        include_prefix + "/pe/kallisto_pe.smk"
    include:
        include_prefix + "/pe/star2_pe.smk"
    include:
        include_prefix + "/pe/reads_feature_count_pe.smk"

