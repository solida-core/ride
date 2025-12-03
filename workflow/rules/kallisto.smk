rule kallisto_index:
    """
    Build kallisto index from transcriptome reference FASTA.
    """
    input:
        fasta=ref_path("transcriptome", "cdna")
    output:
        idx=resolve_results_filepath("kallisto", f"index/{config.get('kallisto').get('index_name')}")
    log:
        resolve_logs_filepath("kallisto","kallisto_index.log")
    threads:
        conservative_cpu_count()
    conda:
        resolve_envs_filepath("kallisto.yaml")
    resources:
        tmpdir=temp_path()
    shell:
        r"""
        mkdir -p $(dirname {output.idx})
        kallisto index \
            -i {output.idx} \
            {input.fasta} \
            > {log} 2>&1
        """

rule kallisto_quant_pe:
    """
    Quantify expression for paired-end libraries with Kallisto.
    """
    input:
        index=rules.kallisto_index.output.idx,
        r1=rules.trim_pe.output.r1,
        r2=rules.trim_pe.output.r2
    output:
        h5=resolve_results_filepath("kallisto", "{sample}/abundance.h5"),
        tsv=resolve_results_filepath("kallisto","{sample}/abundance.tsv"),
        json=resolve_results_filepath("kallisto","{sample}/run_info.json")
    params:
        outdir=resolve_results_filepath("kallisto", "{sample}"),
        boot=config["kallisto"]["boot"]
    log:
        resolve_logs_filepath("kallisto", "{sample}.kallisto.pe.log")
    threads:
        conservative_cpu_count()
    conda:
        resolve_envs_filepath("kallisto.yaml")
    shell:
        r"""
        mkdir -p {params.outdir}
        kallisto quant \
            -i {input.index} \
            -o {params.outdir} \
            -b {params.boot} \
            -t {threads} \
            {input.r1} {input.r2} \
            > {log} 2>&1
        """

rule kallisto_quant_se:
    """
    Quantify expression for single-end libraries with Kallisto.
    """
    input:
        index = rules.kallisto_index.output.idx,
        fastq = rules.trim_se.output.fastq
    output:
        h5=resolve_results_filepath("kallisto", "{sample}/se/abundance.h5"),
        tsv=resolve_results_filepath("kallisto","{sample}/se/abundance.tsv"),
        json=resolve_results_filepath("kallisto","{sample}/se/run_info.json")
    params:
        outdir=resolve_results_filepath("kallisto", "{sample}/se"),
        frag_len = config["kallisto"]["frag_len"],
        frag_sd = config["kallisto"]["frag_sd"],
        boot=config["kallisto"]["boot"]
    log:
        resolve_logs_filepath("kallisto", "{sample}.kallisto.se.log")
    threads:
        conservative_cpu_count()
    conda:
        resolve_envs_filepath("kallisto.yaml")
    shell:
        r"""
        mkdir -p {params.outdir}
        kallisto quant \
            -i {input.index} \
            --single \
            -b {params.boot} \
            -l {params.frag_len} \
            -s {params.frag_sd} \
            -o {params.outdir} \
            -t {threads} \
            {input.fastq} \
            > {log} 2>&1
        """







