import pandas as pd
from snakemake.utils import validate, min_version
import re
##### set minimum snakemake version #####
min_version("5.1.2")


##### load config and sample sheets #####
#configfile: "config.yaml"

## USER FILES ##
samples = pd.read_csv(config["samples"], index_col="sample", sep="\t")
units = pd.read_csv(config["units"], index_col=["unit"], dtype=str, sep="\t")
reheader = pd.read_csv(config["reheader"],index_col="Client", dtype=str, sep="\t")
reheader = reheader[reheader["LIMS"].isin(samples.index.values)]
## ---------- ##


include:
    "rules/functions.py"



# rule fq:
#     input:
#         expand("reads/untrimmed/merged/{sample.sample}-R1.fq.gz",sample=samples.reset_index().itertuples()),
#         expand("reads/untrimmed/merged/{sample.sample}-R2.fq.gz",sample=samples.reset_index().itertuples())

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
        "results/Heatmap_Most_Var.png",
        expand("qc/bbmap_qchist/{sample.sample}-R1.fq.gz.qchist",sample=samples.reset_index().itertuples()),
        expand("delivery/abundance/{Client.Client}.tsv",Client=reheader.reset_index().itertuples()),
        expand("delivery/bams/{Client.Client}.bam",Client=reheader.reset_index().itertuples()),
        expand("delivery/bams/{Client.Client}.bam.bai",Client=reheader.reset_index().itertuples())




include_prefix="rules"
include:
    include_prefix + "/rseqc.smk"
include:
    include_prefix + "/plots.smk"
include:
    include_prefix + "/concatenate_fq.smk"
include:
    include_prefix + "/bbmap.smk"
include:
    include_prefix + "/delivery.smk"
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
    include:
        include_prefix + "/bbduk_se.smk"
else:
    include:
        include_prefix + "/trimming_pe.smk"
    include:
        include_prefix + "/qc_pe.smk"
    include:
        include_prefix + "/kallisto_pe.smk"
    include:
        include_prefix + "/star2_pe.smk"
    include:
        include_prefix + "/reads_feature_count_pe.smk"
    include:
        include_prefix + "/bbduk_pe.smk"
