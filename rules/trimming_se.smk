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
        "reads/untrimmed/merged/{sample}-R1.fq.gz"
    output:
        "reads/trimmed/{sample}-R1_trimmed.fq.gz",
        "reads/trimmed/{sample}-R1.fq.gz_trimming_report.txt"
    params:
        extra=config.get("rules").get("trim_galore_se").get("arguments"),
        outdir="reads/trimmed/"
    log:
        "logs/trim_galore/{sample}.log"
    benchmark:
        "benchmarks/trim_galore/{sample}.txt"
    conda:
        "../envs/trim_galore.yaml"
    threads: (conservative_cpu_count(reserve_cores=2, max_cores=99))/2 if (conservative_cpu_count(reserve_cores=2, max_cores=99)) >2 else 1
    shell:
        "mkdir -p qc/fastqc; "
        "trim_galore "
        "{params.extra} "
        "--cores {threads} "
        "-o {params.outdir} "
        "{input} "
        ">& {log}"




rule post_rename_fastq_se:
    input:
        r1="reads/trimmed/{sample}-R1_trimmed.fq.gz"
    output:
        r1="reads/trimmed/{sample}-R1-trimmed.fq.gz"
    shell:
        "mv {input.r1} {output.r1}"

