
rule fastq_merge_pe_r1:
    input:
        lambda wildcards: get_unit_fastqs_pe(wildcards, read_pair="fq1")
    output:
        fq = resolve_results_filepath(
            "reads",
            "untrimmed/{sample}-R1.fq.gz"
        )
    log:
        resolve_logs_filepath("merge_fq","{sample}-R1.merge.log")
    threads:
        conservative_cpu_count()
    resources:
        tmpdir=temp_path()
    script:
        resolve_scripts_filepath("merge_fq.py")


rule fastq_merge_pe_r2:
    input:
        lambda wildcards: get_unit_fastqs_pe(wildcards, read_pair="fq2")
    output:
        fq = resolve_results_filepath(
            "reads",
            "untrimmed/{sample}-R2.fq.gz"
        )
    log:
        resolve_logs_filepath("merge_fq","{sample}-R2.merge.log")
    threads:
        conservative_cpu_count()
    resources:
        tmpdir=temp_path()
    script:
        resolve_scripts_filepath("merge_fq.py")


rule fastq_merge_se:
    input:
        lambda wildcards: get_unit_fastqs_se(wildcards, read_pair="fq1")
    output:
        fq = resolve_results_filepath(
            "reads",
            "untrimmed/se/{sample}.fq.gz"
        )
    log:
        resolve_logs_filepath("merge_fq","se/{sample}.merge.log")
    threads:
        conservative_cpu_count()
    resources:
        tmpdir=temp_path()
    script:
        resolve_scripts_filepath("merge_fq.py")






