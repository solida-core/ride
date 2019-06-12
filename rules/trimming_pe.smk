def get_fastq(wildcards,samples,read_pair='fq'):
    return samples.loc[wildcards.sample,
                     [read_pair]].dropna()[0]

rule pre_rename_fastq_pe:
    input:
        r1=lambda wildcards: get_fastq(wildcards, samples, read_pair="fq1"),
        r2=lambda wildcards: get_fastq(wildcards, samples, read_pair="fq2")
    output:
        r1="reads/untrimmed/{sample}-R1.fq.gz",
        r2="reads/untrimmed/{sample}-R2.fq.gz"
    shell:
        "ln -s {input.r1} {output.r1} && "
        "ln -s {input.r2} {output.r2}"


rule trim_galore_pe:
    input:
        "reads/untrimmed/{sample}-R1.fq.gz",
        "reads/untrimmed/{sample}-R2.fq.gz"
    output:
        temp("reads/trimmed/{sample}-R1_val_1.fq.gz"),
        "reads/trimmed/{sample}-R1.fq.gz_trimming_report.txt"
        temp("reads/trimmed/{sample}-R2_val_1.fq.gz"),
        "reads/trimmed/{sample}-R2.fq.gz_trimming_report.txt"
    params:
        extra=config.get("rules").get("trim_galore_pe").get("arguments")
    log:
        "logs/trim_galore_pe/{sample}.log"
    benchmark:
        "benchmarks/trim_galore_pe/{sample}.txt"
    wrapper:
        config.get("wrappers").get("trim_galore_pe")

rule post_rename_fastq_pe:
    input:
        r1="reads/trimmed/{sample}-R1_val_1.fq.gz",
        r2="reads/trimmed/{sample}-R2_val_1.fq.gz"
    output:
        r1="reads/trimmed/{sample}-R1-trimmed.fq.gz",
        r2="reads/trimmed/{sample}-R2-trimmed.fq.gz"
    shell:
        "mv {input.r1} {output.r1} && "
        "mv {input.r2} {output.r2}"

