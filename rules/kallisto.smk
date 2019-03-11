rule kallisto_build_index:
    input:
        resolve_single_filepath(*references_abs_path(ref='transcriptome_reference'), config.get("transcriptome_fasta"))
    output: "kallisto/index/transcriptome.kidx.transcripts"
    conda:
        "../envs/kallisto_build_index.yaml"
    shell:
        "kallisto index -i {output} {input}"

rule kallisto_quant:
    input:
        "reads/trimmed/{sample}-R1-trimmed.fq.gz",
        index=rules.kallisto_build_index.output
    output:
        "kallisto/{sample}/abundance.h5",
        "kallisto/{sample}/abundance.tsv",
        "kallisto/{sample}/run_info.json"
    params:
        "kallisto/{sample}"
    conda:
        "../envs/kallisto_quant.yaml"
    threads: pipeline_cpu_count()
    log:
        "logs/kallisto/{sample}.kallisto_quant.log"
    shell:
        "kallisto quant "
        "-i {input.index} "
        "--single "
        "-b 30 "
        "-l 280 "
        "-s 80 "
        "-o {params} "
        "-t {threads} "
        "{input[0]} "
        ">& {log}"




rule sleuth_run:
   input:
       expand("kallisto/{sample}/abundance.h5", sample=config.get('samples')),
       expand("kallisto/{sample}/abundance.tsv", sample=config.get('samples'))
   output:
       dir="kallisto",
       gene_table="kallisto/DEGS/gene_table.txt",
       sleuth_object="kallisto/DEGS/sleuth_object.RData",
       sleuth_table="kallisto/DEGS/sleuth_table.txt"

   conda:
       "../envs/sleuth_run.yaml"
   params:
       classes=lambda wildcards, input: ",".join(config.get('classes').keys()),
       class1=lambda wildcards, input: ",".join(config.get('classes').get('C1')),
       class2=lambda wildcards, input: ",".join(config.get('classes').get('C2')),
       database=config.get("rules").get("sleuth_run").get("database"),
       dataset=config.get("rules").get("sleuth_run").get("dataset"),
       version=config.get("rules").get("sleuth_run").get("version")
   log:
       "logs/sleuth_run.log"

   threads: pipeline_cpu_count()
   script:
       "scripts/sleuth_script.R"