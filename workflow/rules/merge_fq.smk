
rule fastq_merge_pe_r1:
    input:
        lambda wildcards: get_unit_fastqs_pe(wildcards, read_pair="fq1")
    output:
        resolve_results_filepath(
            "reads",
            "untrimmed/{sample}-R1.fq.gz"
        )
    script:
        resolve_scripts_filepath("merge_fq.py")


rule fastq_merge_pe_r2:
    input:
        lambda wildcards: get_unit_fastqs_pe(wildcards, read_pair="fq2")
    output:
        resolve_results_filepath(
            "reads",
            "untrimmed/{sample}-R2.fq.gz"
        )
    script:
        resolve_scripts_filepath("merge_fq.py")


rule fastq_merge_se:
    input:
        lambda wildcards: get_unit_fastqs_se(wildcards, read_pair="fq1")
    output:
        resolve_results_filepath(
            "reads",
            "untrimmed/se/{sample}.fq.gz"
        )
    script:
        resolve_scripts_filepath("merge_fq.py")






