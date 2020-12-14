rule kallisto_build_index:
    input:
        resolve_single_filepath(*references_abs_path(ref='transcriptome_reference'), config.get("transcriptome_fasta"))
    output: "kallisto/index/transcriptome.kidx.transcripts"
    conda:
        "../envs/kallisto.yaml"
    shell:
        "kallisto index -i {output} {input}"





def fastq_input(r1):
    if config.get("read_type")=="se":
        return r1
    else:
        r2=r1.replace("-R1-", "-R2-")
        reads=[r1,r2]
        return " ".join(reads)




rule kallisto_quant:
    input:
        fastq_input("reads/trimmed/{sample}-R1-trimmed.fq.gz"),
        index=rules.kallisto_build_index.output
    output:
        "kallisto/{sample}/abundance.h5",
        "kallisto/{sample}/abundance.tsv",
        "kallisto/{sample}/run_info.json",
        touch("results/tmp/{sample}.ready.for.plots")
    params:
        outdir="kallisto/{sample}",
        params=config.get("rules").get("kallisto").get("arguments" if config.get("read_type")=="pe" else "arguments_se"),
        #read_type=config.get("rules").get("kallisto").get("read_type")
    conda:
        "../envs/kallisto.yaml"
    threads: conservative_cpu_count(reserve_cores=2, max_cores=99)
    log:
        "logs/kallisto/{sample}.kallisto_quant.log"
    shell:
        "kallisto quant "
        "-i {input.index} "
        "-o {params.outdir} "
        "{params.params} "
        "-t {threads} "
        "{input[0]} "
        ">& {log}"

