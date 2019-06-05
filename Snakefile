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
        # Tophat+Cufflinks workflow
#        'assembly/comparison/all.stats',
#        'diffexp/isoform_exp.diff',
        # Kallisto+Sleuth workflow
#        "kallisto/DEGS/gene_table.txt",
        expand("kallisto/{sample.sample}/abundance.tsv", sample=samples.reset_index().itertuples()),
        # Star
        expand("star/{sample.sample}/{sample.sample}.Aligned.sortedByCoord.out.bam", sample=samples.reset_index().itertuples()),
        # Rseqc
#        expand("rseqc/{sample}/{sample}.bam_stat.txt", sample=config.get('samples')),
#        expand("rseqc/{sample}/{sample}.geneBodyCoverage.txt", sample=config.get('samples')),
#        expand("rseqc/{sample}/{sample}.junction.txt", sample=config.get('samples')),
#        expand("rseqc/{sample}/{sample}.junctionSaturation_plot.r", sample=config.get('samples')),
#        expand("rseqc/{sample}/{sample}.GC.xls", sample=config.get('samples')),
##        expand("rseqc/{sample}/{sample}.read_distribution.txt", sample=config.get('samples')),
#        expand("rseqc/{sample}/{sample}.infer_experiment.txt", sample=config.get('samples')),
#        expand("rseqc/{sample}/{sample}.pos.DupRate.xls", sample=config.get('samples')),
#        expand("rseqc/{sample}/{sample}.saturation.pdf", sample=config.get('samples')),
        "qc/multiqc.html"
        # Star+htseq workflow
#        expand("star/{sample}/count/{sample}_counts.cnt", sample=config.get('samples')),


include_prefix="rules"

include:
    include_prefix + "/kallisto.smk"
include:
    include_prefix + "/tophat.smk"
include:
    include_prefix + "/cufflinks.smk"
include:
    include_prefix + "/star2.smk"
include:
    include_prefix + "/rseqc.smk"
include:
    include_prefix + "/qc.smk"
include:
    include_prefix + "/trimming.smk"