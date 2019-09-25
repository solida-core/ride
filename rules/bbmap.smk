
rule bbmap_qchist_r1:
    input:
        "reads/untrimmed/merged/{sample}-R1.fq.gz"
    output:
        "qc/bbmap_qchist/{sample}-R1.fq.gz.qchist"
    conda:
        "../envs/bbmap.yaml"
    log:
        "logs/bbmap/{sample}_qchist.log"
    shell:
        "reformat.sh "
        "in={input} "
        "qchist={output} "
        ">& {log}"

rule bbmap_qchist_r2:
    input:
        "reads/untrimmed/merged/{sample}-R2.fq.gz"
    output:
        "qc/bbmap_qchist/{sample}-R2.fq.gz.qchist"
    conda:
        "../envs/bbmap.yaml"
    log:
        "logs/bbmap/{sample}_qchist.log"
    shell:
        "reformat.sh "
        "in={input} "
        "qchist={output} "
        ">& {log}"