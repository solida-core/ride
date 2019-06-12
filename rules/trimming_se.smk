def get_fastq(wildcards,samples,read_pair='fq'):
    return samples.loc[wildcards.sample,
                     [read_pair]].dropna()[0]

rule pre_rename_fastq_se:
    input:
        r1=lambda wildcards: get_fastq(wildcards, samples, read_pair="fq1")
    output:
        r1="reads/untrimmed/{sample}-R1.fq.gz"
    shell:
        "ln -s {input.r1} {output.r1}"


rule trim_galore_se:
    input:
        "reads/untrimmed/{sample}-R1.fq.gz"
    output:
        temp("reads/trimmed/{sample}-R1_val_1.fq.gz"),
        "reads/trimmed/{sample}-R1.fq.gz_trimming_report.txt"
    params:
        extra=config.get("rules").get("trim_galore_se").get("arguments")
    log:
        "logs/trim_galore/{sample}.log"
    benchmark:
        "benchmarks/trim_galore/{sample}.txt"
    wrapper:
        config.get("wrappers").get("trim_galore_se")

rule post_rename_fastq_pe:
    input:
        r1="reads/trimmed/{sample}-R1_val_1.fq.gz"
    output:
        r1="reads/trimmed/{sample}-R1-trimmed.fq.gz"
    shell:
        "mv {input.r1} {output.r1}"

