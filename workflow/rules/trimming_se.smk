rule pre_rename_fastq_se:
    input:
        r1=lambda wildcards: get_fastq(wildcards, units),
    output:
        r1=resolve_results_filepath("reads", "untrimmed/se/{unit}-R1.fq.gz"),
    log:
        resolve_logs_filepath("bash", "untrimmed/se/{unit}_cp.log"),
    conda:
        resolve_envs_filepath("bash.yaml"),
    resources:
        tmpdir=temp_path(),
    benchmark:
        resolve_benchmarks_filepath('reads',"untrimmed/se/{unit}-R1.txt"),
    threads: conservative_cpu_count(reserve_cores=2,max_cores=99)
    shell:
        "cp {input.r1} {output.r1} "
        ">& {log} "


rule trim_galore_se:
    input:
        r1=rules.pre_rename_fastq_se.output.r1,
    output:
        fq1=resolve_results_filepath(
            "reads", "trimmed/se/{unit}-R1_trimmed.fq.gz"
        ),
        report_fq1=resolve_results_filepath(
            "reads", "trimmed/se/{unit}-R1.fq.gz_trimming_report.txt",
        ),
    params:
        extra=config.get("params").get("trim_galore_se").get("arguments"),
        outdir=lambda w, output: os.path.dirname(output[0]),
        qc_dir=resolve_results_filepath("qc", "fastqc"),
    log:
        resolve_logs_filepath("trim_galore", "/se/{unit}.log"),
    benchmark:
        resolve_benchmarks_filepath("trim_galore", "se/{unit}.txt")
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
        "{input.r1} "
        ">& {log} "

rule post_rename_fastq_se:
    input:
        r1=rules.trim_galore_se.output.fq1,
    output:
        r1=resolve_results_filepath(
            "reads", "/se/trimmed/{unit}-R1-trimmed.fq.gz",
        ),
    log:
        resolve_logs_filepath(
            "bash", "{unit}_mv.log"
        ),
    conda:
        resolve_envs_filepath("bash.yaml"),
    resources:
        tmpdir=temp_path(),
    shell:
        "mv {input} {output.r1}"
        ">& {log} "




