
rule trim_galore_pe:
    input:
        ["reads/untrimmed/merged/{sample}-R1.fq.gz", "reads/untrimmed/merged/{sample}-R2.fq.gz"]
    output:
        temp("reads/trimmed/{sample}-R1_val_1.fq.gz"),
        "reads/trimmed/{sample}-R1.fq.gz_trimming_report.txt",
        temp("reads/trimmed/{sample}-R2_val_2.fq.gz"),
        "reads/trimmed/{sample}-R2.fq.gz_trimming_report.txt"
    params:
        extra=config.get("rules").get("trim_galore_pe").get("arguments"),
        outdir="reads/trimmed/"
    log:
        "logs/trim_galore/{sample}.log"
    benchmark:
        "benchmarks/trim_galore/{sample}.txt"
    conda:
        "../envs/trim_galore.yaml"
    threads: (conservative_cpu_count(reserve_cores=2, max_cores=99))/4 if (conservative_cpu_count(reserve_cores=2, max_cores=99)) >4 else 1
    shell:
        "mkdir -p qc/fastqc; "
        "trim_galore "
        "{params.extra} "
        "--cores {threads} "
        "-o {params.outdir} "
        "{input} "
        ">& {log}"





rule post_rename_fastq_pe:
    input:
        r1="reads/trimmed/{sample}-R1_val_1.fq.gz",
        r2="reads/trimmed/{sample}-R2_val_2.fq.gz"
    output:
        r1="reads/trimmed/{sample}-R1-trimmed.fq.gz",
        r2="reads/trimmed/{sample}-R2-trimmed.fq.gz"
    shell:
        "mv {input.r1} {output.r1} && "
        "mv {input.r2} {output.r2}"

