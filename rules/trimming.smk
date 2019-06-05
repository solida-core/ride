
rule pre_rename_fastq_pe:
    input:
        r1=lambda wildcards: config["samples"][wildcards.sample]
    output:
        r1="reads/untrimmed/{sample}-R1.fq.gz"
    shell:
        "ln -s {input.r1} {output.r1}"


rule trim_galore_pe:
    input:
        "reads/{sample}-R1.fq.gz"
    output:
        temp("reads/trimmed/{sample}-R1_val_1.fq.gz"),
        "reads/trimmed/{sample}-R1.fq.gz_trimming_report.txt"
    params:
        extra=config.get("rules").get("trim_galore_pe").get("arguments")
    log:
        "logs/trim_galore/{sample}.log"
    benchmark:
        "benchmarks/trim_galore/{sample}.txt"
    wrapper:
#        "0.27.0/bio/trim_galore/pe"
        config.get("wrappers").get("trim_galore")

rule post_rename_fastq_pe:
    input:
        r1="reads/trimmed/{sample}-R1_val_1.fq.gz"
    output:
        r1="reads/trimmed/{sample}-R1-trimmed.fq.gz"
    shell:
        "mv {input.r1} {output.r1}"

