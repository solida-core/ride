rule fastqc_pe:
    """
    Run FastQC on paired-end trimmed reads.
    """
    input:
        r1 = rules.trim_pe.output.r1,
        r2 = rules.trim_pe.output.r2
    output:
        html_r1 = resolve_results_filepath(
            "qc",
            "fastqc/{sample}-R1.trimmed_fastqc.html"
        ),
        zip_r1  = resolve_results_filepath(
            "qc",
            "fastqc/{sample}-R1.trimmed_fastqc.zip"
        ),
        html_r2 = resolve_results_filepath(
            "qc",
            "fastqc/{sample}-R2.trimmed_fastqc.html"
        ),
        zip_r2  = resolve_results_filepath(
            "qc",
            "fastqc/{sample}-R2.trimmed_fastqc.zip"
        )
    params:
        outdir = resolve_results_filepath(
            "qc",
            "fastqc"
        )
    log:
        resolve_logs_filepath(
            "fastqc",
            "{sample}.fastqc.pe.log"
        )
    benchmark:
        resolve_benchmarks_filepath(
            "fastqc",
            "{sample}.fastqc.pe.txt"
        )
    threads:
        conservative_cpu_count()
    conda:
        resolve_envs_filepath("quality_check.yaml")
    resources:
        tmpdir=temp_path()
    shell:
        # "mkdir -p {params.outdir} ; "
        "fastqc "
        "--threads {threads} "
        "--outdir {params.outdir} "
        "{input.r1} {input.r2} "
        ">& {log} "


rule fastqc_se:
    """
    Run FastQC on single-end trimmed reads.
    """
    input:
        fq = rules.trim_se.output.fastq
    output:
        html = resolve_results_filepath(
            "qc",
            "fastqc/se/{sample}.trimmed_fastqc.html"
        ),
        zip  = resolve_results_filepath(
            "qc",
            "fastqc/se/{sample}.trimmed_fastqc.zip"
        )
    params:
        outdir=resolve_results_filepath(
            "qc",
            "fastqc/se"
        )
    log:
        resolve_logs_filepath(
            "fastqc",
            "{sample}.fastqc.se.log"
        )
    benchmark:
        resolve_benchmarks_filepath(
            "fastqc",
            "{sample}.fastqc.se.txt"
        )
    threads:
        conservative_cpu_count()
    conda:
        resolve_envs_filepath("quality_check.yaml")
    resources:
        tmpdir=temp_path()
    shell:
        "fastqc "
        "--threads {threads} "
        "--outdir {params.outdir} "
        "{input.fq} "
        ">& {log} "


rule multiqc:
    """
    Run MultiQC on all QC outputs (FastQC + Fastp + Salmon/Kallisto).
    """
    input:
        # FastQC
        expand(
            resolve_results_filepath(
                "qc",
                "fastqc/{sample}-R1.trimmed_fastqc.zip"),
            sample=SAMPLES_PE
        ),
        expand(
            resolve_results_filepath(
                "qc",
                "fastqc/{sample}-R2.trimmed_fastqc.zip"),
            sample=SAMPLES_PE
        ),
        expand(
            resolve_results_filepath(
                "qc",
                "fastqc/se/{sample}.trimmed_fastqc.zip"),
            sample=SAMPLES_SE
        ),

        # Trimmimg - Fastp
        expand(
            resolve_results_filepath(
                "qc",
                "trimming/{sample}.fastp.pe.json"),
            sample=SAMPLES_PE
        ),
        expand(
            resolve_results_filepath(
                "qc",
                "trimming/se/{sample}.fastp.se.json"),
            sample=SAMPLES_SE
        ),

        # Kallisto
        expand(
            resolve_logs_filepath(
                "kallisto",
                "{sample}.kallisto.pe.log"),
            sample=SAMPLES_PE
        ),
        expand(
            resolve_logs_filepath(
                "kallisto",
                "{sample}.kallisto.se.log"),
            sample=SAMPLES_SE
        ),

        # Star
        expand(
            resolve_results_filepath(
                "star",
                "{sample}/{sample}.pe.Log.final.out"),
            sample=SAMPLES_PE
        ),
        expand(
            resolve_results_filepath(
                "star",
                "{sample}/se/{sample}.se.Log.final.out"),
            sample=SAMPLES_SE
        )
    output:
        html = resolve_results_filepath(
            "qc",
            "multiqc/multiqc_report.html"
        ),
        data = directory(
            resolve_results_filepath(
                "qc",
                "multiqc/multiqc_data"
            )
        )
    params:
        outdir = resolve_results_filepath(
            "qc",
            "multiqc"
        ),
        fastqc = resolve_results_filepath(
            "qc",
            "fastqc"
        ),
        trimming =  resolve_results_filepath(
            "qc",
            "trimming"
        ),
        kallisto = resolve_logs_filepath(
            "kallisto",
            ""
        ),
        star = resolve_logs_filepath(
            "star",
            ""
        ),
    log:
        resolve_logs_filepath(
            "multiqc",
            "multiqc.log"
        )
    threads:
        conservative_cpu_count()
    conda:
        resolve_envs_filepath("quality_check.yaml")
    shell:
        "multiqc "
        "{params.star} "
        "{params.kallisto} "
        "{params.trimming} "
        "{params.fastqc} "
        "-o {params.outdir} "
        ">& {log} "

