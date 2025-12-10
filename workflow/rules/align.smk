rule star_index:
    """
    Build STAR genome index from reference FASTA and GTF.
    """
    input:
        fasta = ref_path(
            "genome",
            "fasta"
        ),
        gtf   = ref_path(
            "genome",
            "gtf")
    output:
        indexdir = directory(
            resolve_results_filepath(
                "star",
                config["star"]["index_dir"]
            )
        ),
        length = resolve_results_filepath(
            "star",
            f"{config['star']['index_dir']}/chrLength.txt"
        )

    params:
        indexdir = resolve_results_filepath(
            "star",
            config["star"]["index_dir"]
        ),
        sjdb_overhang = config["star"]["sjdb_overhang"]
    log:
        resolve_logs_filepath(
            "star",
            "star_index.log"
        )
    benchmark:
        resolve_benchmarks_filepath(
            "star",
            "star_index.txt"
        )
    threads:
        conservative_cpu_count()
    conda:
        resolve_envs_filepath("align.yaml")
    resources:
        tmpdir = temp_path()
    shell:
        "STAR "
        "--runMode genomeGenerate "
        "--genomeDir {params.indexdir} "
        "--genomeFastaFiles {input.fasta} "
        "--sjdbGTFfile {input.gtf} "
        "--sjdbOverhang {params.sjdb_overhang} "
        "--runThreadN {threads} "
        ">& {log} "

rule star_align_pe:
    """
    Align paired-end trimmed reads with STAR and produce sorted BAM.
    """
    input:
        idx = rules.star_index.output.indexdir,
        r1  = rules.trim_pe.output.r1,
        r2  = rules.trim_pe.output.r2,
        length = rules.star_index.output.length
    output:
        bam = resolve_results_filepath(
            "star",
            "{sample}/{sample}.bam"
        )
    params:
        outprefix = resolve_results_filepath(
            "star",
            "{sample}/{sample}.pe."
        ),
        extra     = config["star"].get("extra_params", "")
    log:
        resolve_logs_filepath(
            "star",
            "{sample}.star.pe.log"
        )
    benchmark:
        resolve_benchmarks_filepath(
            "star",
            "{sample}.star.pe.txt"
        )
    threads:
        conservative_cpu_count()
    conda:
        resolve_envs_filepath("align.yaml")
    resources:
        tmpdir = temp_path()
    shell:
        "STAR "
        "--genomeDir {input.idx} "
        "--readFilesIn {input.r1} {input.r2} "
        "--readFilesCommand zcat "
        "--runThreadN {threads} "
        "--outFileNamePrefix {params.outprefix} "
        "--outSAMtype BAM SortedByCoordinate "
        "--outSAMunmapped Within "
        "{params.extra} "
        ">& {log} ; "
        "mv  {params.outprefix}Aligned.sortedByCoord.out.bam {output.bam}"

rule star_align_se:
    """
    Align single-end trimmed reads with STAR and produce sorted BAM.
    """
    input:
        idx  = rules.star_index.output.indexdir,
        fastq = rules.trim_se.output.fastq
    output:
        bam = resolve_results_filepath(
            "star",
            "{sample}/se/{sample}.bam"
        )
    params:
        outprefix = resolve_results_filepath(
            "star",
            "{sample}/se/{sample}.se."
        ),
        extra     = config["star"].get("extra_params", "")
    log:
        resolve_logs_filepath(
            "star",
            "{sample}.star.se.log")
    benchmark:
        resolve_benchmarks_filepath(
            "star",
            "{sample}.star.se.txt"
        )
    threads:
        conservative_cpu_count()
    conda:
        resolve_envs_filepath("align.yaml")
    resources:
        tmpdir = temp_path()
    shell:
        "STAR "
        "--genomeDir {input.idx} "
        "--readFilesIn {input.fastq} "
        "--readFilesCommand zcat "
        "--runThreadN {threads} "
        "--outFileNamePrefix {params.outprefix} "
        "--outSAMtype BAM SortedByCoordinate "
        "--outSAMunmapped Within "
        "--quantMode GeneCounts "
        "{params.extra} "
        ">& {log} ; "
        "mv  {params.outprefix}Aligned.sortedByCoord.out.bam {output.bam}"

rule index_bam_pe:
    input:
        bam = rules.star_align_pe.output.bam
    output:
        bai = resolve_results_filepath(
            "star",
            "{sample}/{sample}.bam.bai"
        )
    conda:
        resolve_envs_filepath("align.yaml")
    log:
        resolve_logs_filepath(
            "star",
            "{sample}.index.pe.log"
        )
    benchmark:
        resolve_benchmarks_filepath(
            "star",
            "{sample}.index.pe.txt"
        )
    threads:
        conservative_cpu_count()
    shell:
        "samtools index "
        "--threads {threads} "
        "{input.bam} {output.bai} "
        ">& {log} "


rule index_bam_se:
    input:
        bam = rules.star_align_se.output.bam
    output:
        bai = resolve_results_filepath(
            "star",
            "{sample}/se/{sample}.bam.bai"
        )
    conda:
        resolve_envs_filepath("align.yaml")
    log:
        resolve_logs_filepath(
            "star",
            "{sample}.index.se.log"
        )
    benchmark:
        resolve_benchmarks_filepath(
            "star",
            "{sample}.index.se.txt"
        )
    threads:
        conservative_cpu_count()
    shell:
        "samtools index "
        "--threads {threads} "
        "{input.bam} {output.bai} "
        ">& {log} "


