
rule multiqc:
    input:
multiqc -m kallisto -m star -m rseqc ...

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