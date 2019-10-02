

rule featureCounts_run:
    input:
        "reads/trimmed/{sample}-R1-trimmed.fq.gz",
        "reads/trimmed/{sample}-R2-trimmed.fq.gz",
        bam="star/{sample}/{sample}.Aligned.sortedByCoord.out.bam"
    output:
        "star/{sample}/count/{sample}_featurecounts.cnt"
    conda:
        "../envs/featureCounts.yaml"
    params:
        cmd="featureCounts",
        gtf=resolve_single_filepath(*references_abs_path(ref='references'),
                                    config.get("genes_gtf")),
        gtf_feature_type=config.get("rules").get("featureCounts_run").get("gtf_feature_type"),
        gtf_attribute_type=config.get("rules").get("featureCounts_run").get("gtf_attribute_type"),
    threads: pipeline_cpu_count()
    script:
        "../scripts/featureCounts_script.py"


rule HTSeq_run:
    input:
        "reads/trimmed/{sample}-R1-trimmed.fq.gz",
        "reads/trimmed/{sample}-R2-trimmed.fq.gz",
        bam="star/{sample}/{sample}.Aligned.sortedByCoord.out.bam"
    output:
        "star/{sample}/count/{sample}_HTSeqcounts.cnt"
    conda:
        "envs/htseq.yaml"
    log:
        "star/{sample}/log/{sample}_htseq_count.log"
    params:
         gtf=resolve_single_filepath(*references_abs_path(ref='references'),
                                    config.get("genes_gtf")),
         strand=config['strand'],
    shell:
         "htseq-count "
         "-m intersection-nonempty "
         "--stranded={params.strand} "
         "--idattr gene_id "
         "-r pos "
         "-f bam "
         "{input.bam} {params.gtf} > {output} 2> {log}"
