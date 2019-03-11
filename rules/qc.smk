
rule multiqc:
    input:
         expand("qc/untrimmed_{sample}.html", sample=config.get('samples')),
         expand("qc/trimmed_{sample}.html", sample=config.get('samples')),
         expand("reads/trimmed/{sample}-R1.fq.gz_trimming_report.txt", sample=config.get('samples')),
         expand("rseqc/{sample}/{sample}.bam_stat.txt", sample=config.get('samples')),
         expand("rseqc/{sample}/{sample}.geneBodyCoverage.txt", sample=config.get('samples')),
         expand("rseqc/{sample}/{sample}.junction.txt", sample=config.get('samples')),
         expand("rseqc/{sample}/{sample}.junctionSaturation_plot.r", sample=config.get('samples')),
         expand("rseqc/{sample}/{sample}.GC.xls", sample=config.get('samples')),
         expand("rseqc/{sample}/{sample}.read_distribution.txt", sample=config.get('samples')),
         expand("rseqc/{sample}/{sample}.infer_experiment.txt", sample=config.get('samples')),
         expand("rseqc/{sample}/{sample}.pos.DupRate.xls", sample=config.get('samples')),
         expand("star/{sample}/{sample}.Log.final.out", sample=config.get('samples')),
         expand("logs/kallisto/{sample}.kallisto_quant.log", sample=config.get('samples'))
    output:
        "qc/multiqc.html"
    params:
        params=config.get("rules").get("multiqc").get("arguments"),
        outdir="qc",
        outname="multiqc.html"
    conda:
        "../envs/multiqc.yaml"
    log:
        "logs/multiqc/multiqc.log"
    shell:
        "multiqc "
        "{input} "
        "{params.params} "
        "-o {params.outdir} "
        "-n {params.outname} "
        ">& {log}"

rule fastqc:
    input:
       "reads/{sample}-R1.fq.gz"
    output:
        html="qc/untrimmed_{sample}.html",
        zip="qc/untrimmed_{sample}_fastqc.zip"
    log:
        "logs/fastqc/{sample}.log"
    params: ""
    wrapper:
        config.get("wrappers").get("fastqc")

rule fastqc_trimmed:
    input:
       "reads/trimmed/{sample}-R1-trimmed.fq.gz"
    output:
        html="qc/trimmed_{sample}.html",
        zip="qc/trimmed_{sample}_fastqc.zip"
    log:
        "logs/fastqc/{sample}.log"
    params: ""
    wrapper:
        config.get("wrappers").get("fastqc")