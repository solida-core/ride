
rule all:
    input:
        # Tophat+Cufflinks workflow
#        'assembly/comparison/all.stats',
#        'diffexp/isoform_exp.diff',
        # Kallisto+Sleuth workflow
#        "kallisto/DEGS/gene_table.txt",
        expand("kallisto/{sample}/abundance.tsv", sample=config.get('samples')),
        # Star
        expand("star/{sample}/{sample}.Aligned.sortedByCoord.out.bam", sample=config.get('samples')),
        # Rseqc
        expand("rseqc/{sample}/{sample}.bam_stat.txt", sample=config.get('samples')),
        expand("rseqc/{sample}/{sample}.geneBodyCoverage.txt", sample=config.get('samples')),
        expand("rseqc/{sample}/{sample}.junction.txt", sample=config.get('samples')),
        expand("rseqc/{sample}/{sample}.junctionSaturation_plot.r", sample=config.get('samples')),
        expand("rseqc/{sample}/{sample}.GC.xls", sample=config.get('samples')),
        expand("rseqc/{sample}/{sample}.read_distribution.txt", sample=config.get('samples')),
        expand("rseqc/{sample}/{sample}.infer_experiment.txt", sample=config.get('samples')),
        expand("rseqc/{sample}/{sample}.pos.DupRate.xls", sample=config.get('samples')),
        expand("rseqc/{sample}/{sample}.saturation.pdf", sample=config.get('samples')),
        "qc/multiqc.html"
        # Star+htseq workflow
#        expand("star/{sample}/count/{sample}_counts.cnt", sample=config.get('samples')),


include_prefix="rules"

include:
    include_prefix + "/functions.py"
include:
    include_prefix + "/kallisto.smk"
include:
    include_prefix + "/tophat.smk"
include:
    include_prefix + "/cufflinks.smk"
include:
    include_prefix + "/notify.smk"
include:
    include_prefix + "/star2.smk"
include:
    include_prefix + "/rseqc.smk"