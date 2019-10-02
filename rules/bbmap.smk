
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

rule qc_hist_2tsv:
    input:
        expand("qc/bbmap_qchist/{sample.sample}-R1.fq.gz.qchist",sample=samples.reset_index().itertuples())
    output:
        table="qc/bbmap_qchist_summary.tsv"
    params:
        current_dir=get_cwd(),
        input_path="qc/bbmap_qchist",
        out_dir="qc/"
    conda:
        "../envs/rplots.yaml"
    script:
        "scripts/qchist_to_table.R"