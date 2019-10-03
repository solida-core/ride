rule kallisto_build_index:
    input:
        resolve_single_filepath(*references_abs_path(ref='transcriptome_reference'), config.get("transcriptome_fasta"))
    output: "kallisto/index/transcriptome.kidx.transcripts"
    conda:
        "../envs/kallisto.yaml"
    shell:
        "kallisto index -i {output} {input}"



rule kallisto_quant:
    input:
        "reads/trimmed/{sample}-R1-trimmed.fq.gz",
        "reads/trimmed/{sample}-R2-trimmed.fq.gz",
        index=rules.kallisto_build_index.output
    output:
        "kallisto/{sample}/abundance.h5",
        "kallisto/{sample}/abundance.tsv",
        "kallisto/{sample}/run_info.json",
        touch("results/tmp/{sample}.ready.for.plots")
    params:
        outdir="kallisto/{sample}",
        params=config.get("rules").get("kallisto").get("arguments_pe" if config.get("read_type")=="pe" else "arguments_se"),
        #read_type=config.get("rules").get("kallisto").get("read_type")
    conda:
        "../envs/kallisto.yaml"
    threads: pipeline_cpu_count()
    log:
        "logs/kallisto/{sample}.kallisto_quant.log"
    shell:
        "kallisto quant "
        "-i {input.index} "
        "-o {params.outdir} "
        "{params.params} "
        "-t {threads} "
        "{input[0]} {input[1]} "
        ">& {log}"




rule sleuth_run:
   input:
       expand("kallisto/{sample.sample}/abundance.h5", sample=samples.reset_index().itertuples()),
       expand("kallisto/{sample.sample}/abundance.tsv", sample=samples.reset_index().itertuples())
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
       "../scripts/sleuth_script.R"
