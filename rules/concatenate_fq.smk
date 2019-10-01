def get_unit_fastqs(wildcards, samples, label='units',read_pair='fq'):
    for unit_set in samples.loc[wildcards.sample,[label]]:
        print(wildcards.sample)
    return [units.loc[x,[read_pair]].dropna()[0] for x in unit_set.split(',')]

rule fastq_merge_r1:
    input:
        lambda wildcards: get_unit_fastqs(wildcards, samples, read_pair='fq1')
    output:
        "reads/untrimmed/merged/{sample}-R1.fq.gz"
    script:
        "scripts/merge_units.py"


rule fastq_merge_r2:
    input:
        lambda wildcards: get_unit_fastqs(wildcards, samples, read_pair='fq2')
    output:
        "reads/untrimmed/merged/{sample}-R2.fq.gz"
    script:
        "scripts/merge_units.py"




