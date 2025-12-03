
rule fastq_merge_r1:
    input:
        lambda wildcards: get_unit_fastqs(wildcards, samples, read_pair='fq1')
    output:
        "reads/untrimmed/merged/{sample}-R1.fq.gz"
    script:
        resolve_scripts_filepath("merge_units.py")



rule fastq_merge_r2:
    input:
        lambda wildcards: get_unit_fastqs(wildcards, samples, read_pair='fq2')
    output:
        resolve_results_filepath(
            "reads", "untrimmed/merged/{sample}-R2.fq.gz"
        ),
    script:
        resolve_scripts_filepath("merge_units.py")




