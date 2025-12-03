
rule trim_galore_pe:
    input:
        rules.fastq_merge_r1.output,
        rules.fastq_merge_r2.output,
    output:
        fq1=resolve_results_filepath(
            "reads", "trimmed/{unit}-R1_val_1.fq.gz",
        ),
        report_r1=resolve_results_filepath(
            "reads", "trimmed/{unit}-R1.fq.gz_trimming_report.txt",
        ),
        fq2=resolve_results_filepath(
            "reads", "trimmed/{unit}-R2_val_2.fq.gz",
        ),
        report_r2=resolve_results_filepath(
            "reads", "trimmed/{unit}-R2.fq.gz_trimming_report.txt",
        ),
    params:
        extra=config.get("params").get("trim_galore_pe").get("arguments"),
        outdir=lambda w, output: os.path.dirname(output[0]),
        qc_dir=resolve_results_filepath(
            "qc", "fastqc"
        ),
    log:
        resolve_logs_filepath("trim_galore","{sample}.log"),
    benchmark:
        resolve_benchmarks_filepath("trim_galore","{sample}.txt")
    conda:
        resolve_envs_filepath("trim_galore.yaml"),
    resources:
        tmpdir=temp_path(),
    threads: conservative_cpu_count(reserve_cores=2,max_cores=99)
    shell:
        "mkdir -p {params.qc_dir}; "
        "trim_galore "
        "{params.extra} "
        "--cores {threads} "
        "-o {params.outdir} "
        "{input} "
        ">& {log}"

rule post_rename_fastq_pe:
    input:
        r1=rules.trim_galore_pe.output.fq1,
        r2=rules.trim_galore_pe.output.fq2,
    output:
        r1=resolve_results_filepath(
            "reads", "trimmed/{sample}-R1-trimmed.fq.gz",
        ),
        r2=resolve_results_filepath(
            "reads", "trimmed/{sample}-R2-trimmed.fq.gz"
        ),
    log:
        resolve_logs_filepath("bash","pe/{unit}_mv.log"),
    conda:
        resolve_envs_filepath("bash.yaml"),
    resources:
        tmpdir=temp_path(),
    threads: conservative_cpu_count(reserve_cores=2,max_cores=99)
    shell:
        "mv {input.r1} {output.r1} && "
        "mv {input.r2} {output.r2} "
        ">& {log} "
