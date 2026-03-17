rule trim_pe:
    """
      Perform paired-end read trimming with fastp.
      Generates:
        - trimmed R1 and R2 fastq.gz
        - HTML QC report
        - JSON summary
    """
    input:
        r1 = rules.fastq_merge_pe_r1.output,
        r2 = rules.fastq_merge_pe_r2.output
    output:
        r1 = resolve_results_filepath(
            "reads",
            "trimmed/{sample}-R1.trimmed.fq.gz"
        ),
        r2 = resolve_results_filepath(
            "reads",
            "trimmed/{sample}-R2.trimmed.fq.gz"
        ),
        html = resolve_results_filepath(
            "qc",
            "trimming/{sample}.fastp.pe.html"
        ),
        json = resolve_results_filepath(
            "qc",
            "trimming/{sample}.fastp.pe.json"
        )
    params:
        quality = config["trimming"]["quality"],
        min_length = config["trimming"]["min_length"],
        trim_poly_g = "--trim_poly_g" if config["trimming"]["trim_poly_g"] else "",
        trim_poly_x = "--trim_poly_x" if config["trimming"]["trim_poly_x"] else "",
        detect_adapter = "--detect_adapter_for_pe" if config["trimming"]["detect_adapter"] else ""
    log:
        resolve_logs_filepath(
            "trimming",
            "{sample}.fastp.pe.log"
        )
    benchmark:
        resolve_benchmarks_filepath(
            "trimming",
            "{sample}.fastp.pe.txt"
        )
    threads:
        conservative_cpu_count()
    conda:
        resolve_envs_filepath("fastp.yaml")
    resources:
        tmpdir=temp_path()
    shell:
        "fastp "
        "-i {input.r1} "
        "-I {input.r2} "
        "-o {output.r1} "
        "-O {output.r2} "
        "--thread {threads} "
        "--qualified_quality_phred {params.quality} "
        "--length_required {params.min_length} "
        "{params.trim_poly_g} "
        "{params.trim_poly_x} "
        "{params.detect_adapter} "
        "--html {output.html} "
        "--json {output.json} "
        ">& {log} "

rule trim_se:
    """
      Perform single-end read trimming with fastp.
      Generates:
        - trimmed R1 and R2 fastq.gz
        - HTML QC report
        - JSON summary
    """
    input:
        rules.fastq_merge_se.output
    output:
        fastq = resolve_results_filepath(
            "reads",
            "trimmed/se/{sample}.trimmed.fq.gz"
        ),
        html = resolve_results_filepath(
            "qc",
            "trimming/se/{sample}.fastp.se.html"
        ),
        json = resolve_results_filepath(
            "qc",
            "trimming/se/{sample}.fastp.se.json"
        )
    params:
        quality = config["trimming"]["quality"],
        min_length = config["trimming"]["min_length"],
        trim_poly_g = "--trim_poly_g" if config["trimming"]["trim_poly_g"] else "",
        trim_poly_x = "--trim_poly_x" if config["trimming"]["trim_poly_x"] else "",
    log:
        resolve_logs_filepath(
            "trimming",
            "{sample}.fastp.se.log"
        )
    benchmark:
        resolve_benchmarks_filepath(
            "trimming",
            "{sample}.fastp.se.txt"
        )
    threads:
        conservative_cpu_count()
    conda:
        resolve_envs_filepath("fastp.yaml")
    resources:
        tmpdir=temp_path()
    shell:
        "fastp "
        "-i {input} "
        "-o {output.fastq} "
        "--thread {threads} "
        "--qualified_quality_phred {params.quality} "
        "--length_required {params.min_length} "
        "{params.trim_poly_g} "
        "{params.trim_poly_x} "
        "--html {output.html} "
        "--json {output.json} "
        ">& {log} "