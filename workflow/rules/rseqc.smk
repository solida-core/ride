rule rseqc_infer_experiment_pe:
    """
    RSeQC: infer library strandedness for paired-end BAM.
    """
    input:
        bam  = rules.star_align_pe.output.bam,
        bai  = rules.index_bam_pe.output.bai
    output:
        txt  = resolve_results_filepath(
            "rseqc",
            "{sample}/{sample}.pe.infer_experiment.txt"
        )
    params:
        bed = ref_path("rseqc","refseq")
    benchmark:
        resolve_benchmarks_filepath(
            "rseqc",
            "{sample}.rseqc.infer_experiment.pe.txt"
        )
    threads:
        conservative_cpu_count()
    conda:
        resolve_envs_filepath("rseqc.yaml")
    shell:
        "infer_experiment.py "
        "-r {params.bed} "
        "-i {input.bam} "
        ">& {output.txt} "

rule rseqc_infer_experiment_se:
    """
    RSeQC: infer library strandedness for single-end BAM.
    """
    input:
        bam  = rules.star_align_se.output.bam,
        bai  = rules.index_bam_se.output.bai
    output:
        txt  = resolve_results_filepath(
            "rseqc",
            "{sample}/se/{sample}.se.infer_experiment.txt"
        )
    params:
        bed = ref_path("rseqc","bed")
    benchmark:
        resolve_benchmarks_filepath(
            "rseqc",
            "{sample}.rseqc.infer_experiment.pe.txt"
        )
    threads:
        conservative_cpu_count()
    conda:
        resolve_envs_filepath("rseqc.yaml")
    shell:
        "infer_experiment.py "
        "-r {params.bed} "
        "-i {input.bam} "
        ">& {output.txt} "

rule rseqc_read_distribution_pe:
    """
    RSeQC: read distribution over genomic features for paired-end BAM.
    """
    input:
        bam  = rules.star_align_pe.output.bam,
        bai  = rules.index_bam_pe.output.bai,
    output:
        txt  = resolve_results_filepath(
            "rseqc",
            "{sample}/{sample}.pe.read_distribution.txt"
        )
    params:
        bed = ref_path("rseqc","bed")
    benchmark:
        resolve_benchmarks_filepath(
            "rseqc",
            "{sample}.rseqc.read_distribution.pe.txt"
        )
    threads:
        conservative_cpu_count()
    conda:
        resolve_envs_filepath("rseqc.yaml")
    shell:
        "read_distribution.py "
        "-r {params.bed} "
        "-i {input.bam} "
        ">& {output.txt} "

rule rseqc_read_distribution_se:
    """
    RSeQC: read distribution over genomic features for single-end BAM.
    """
    input:
        bam  = rules.star_align_se.output.bam,
        bai  = rules.index_bam_se.output.bai,
    output:
        txt  = resolve_results_filepath(
            "rseqc",
            "{sample}/se/{sample}.se.read_distribution.txt"
        )
    params:
        bed = ref_path("rseqc","bed")
    benchmark:
        resolve_benchmarks_filepath(
            "rseqc",
            "{sample}.rseqc.read_distribution.se.txt"
        )
    threads:
        conservative_cpu_count()
    conda:
        resolve_envs_filepath("rseqc.yaml")
    shell:
        "read_distribution.py "
        "-r {params.bed} "
        "-i {input.bam} "
        ">& {output.txt} "

rule rseqc_geneBody_coverage_pe:
    """
    RSeQC: gene body coverage (5'–3' bias) for paired-end BAM.
    """
    input:
        bam  = rules.star_align_pe.output.bam,
        bai  = rules.index_bam_pe.output.bai,
    output:
        txt = resolve_results_filepath(
            "rseqc",
            "{sample}/{sample}.pe.geneBodyCoverage.txt"
        )
    params:
        bed  = ref_path("rseqc", "bed"),
        out_prefix = resolve_results_filepath(
            "rseqc",
            "{sample}/{sample}.pe"
        )
    benchmark:
        resolve_benchmarks_filepath(
            "rseqc",
            "{sample}.rseqc.geneBody_coverage.pe.txt"
        )
    threads:
        conservative_cpu_count()
    conda:
        resolve_envs_filepath("rseqc.yaml")
    shell:
        "geneBody_coverage.py "
        "-r {params.bed} "
        "-i {input.bam} "
        "-o {params.out_prefix} "

rule rseqc_geneBody_coverage_se:
    """
    RSeQC: gene body coverage (5'–3' bias) for single-end BAM.
    """
    input:
        bam = rules.star_align_se.output.bam,
        bai = rules.index_bam_se.output.bai,
    output:
        txt = resolve_results_filepath(
            "rseqc",
            "{sample}/se/{sample}.se.geneBodyCoverage.txt"
        )
    params:
        bed = ref_path("rseqc","bed"),
        out_prefix=resolve_results_filepath(
            "rseqc",
            "{sample}/se/{sample}.se"
        )
    benchmark:
        resolve_benchmarks_filepath(
            "rseqc",
            "{sample}.rseqc.geneBody_coverage.se.txt"
        )
    threads:
        conservative_cpu_count()
    conda:
        resolve_envs_filepath("rseqc.yaml")
    shell:
        "geneBody_coverage.py "
        "-r {params.bed} "
        "-i {input.bam} "
        "-o {params.out_prefix} "

rule rseqc_bam_stat_pe:
    """
    RSeQC: general BAM statistics for paired-end BAM.
    """
    input:
        bam = rules.star_align_pe.output.bam,
        bai = rules.index_bam_pe.output.bai
    output:
        txt = resolve_results_filepath(
            "rseqc",
            "{sample}/{sample}.pe.bam_stat.txt"
        )
    benchmark:
        resolve_benchmarks_filepath(
            "rseqc",
            "{sample}.rseqc.bam_stat.pe.txt"
        )
    threads:
        conservative_cpu_count()
    conda:
        resolve_envs_filepath("rseqc.yaml")
    shell:
        "bam_stat.py "
        "-i {input.bam} "
        ">& {output.txt} "

rule rseqc_bam_stat_se:
    """
    RSeQC: general BAM statistics for single-end BAM.
    """
    input:
        bam = rules.star_align_se.output.bam,
        bai = rules.index_bam_se.output.bai
    output:
        txt = resolve_results_filepath(
            "rseqc",
            "{sample}/se/{sample}.se.bam_stat.txt"
        )
    benchmark:
        resolve_benchmarks_filepath(
            "rseqc",
            "{sample}.rseqc.bam_stat.se.txt"
        )
    threads:
        conservative_cpu_count()
    conda:
        resolve_envs_filepath("rseqc.yaml")
    shell:
        "bam_stat.py "
        "-i {input.bam} "
        ">& {output.txt} "

rule rseqc_junction_annotation_pe:
    """
    RSeQC: splice junction annotation for paired-end BAM.
    """
    input:
        bam = rules.star_align_pe.output.bam,
        bai = rules.index_bam_pe.output.bai,
    output:
        txt = resolve_results_filepath(
            "rseqc",
            "{sample}/{sample}.pe.junction_annotation.txt"
        )
    params:
        bed = ref_path("rseqc", "bed"),
        out_prefix = resolve_results_filepath(
            "rseqc",
            "{sample}/{sample}.pe."
        )
    benchmark:
        resolve_benchmarks_filepath(
            "rseqc",
            "{sample}.rseqc.junction_annotation.pe.txt"
        )
    threads:
        conservative_cpu_count()
    conda:
        resolve_envs_filepath("rseqc.yaml")
    shell:
        "junction_annotation.py "
        "-r {params.bed} "
        "-i {input.bam} "
        "-o {params.out_prefix} "
        ">& {output.txt} "

rule rseqc_junction_annotation_se:
    """
    RSeQC: splice junction annotation for single-end BAM.
    """
    input:
        bam = rules.star_align_se.output.bam,
        bai = rules.index_bam_se.output.bai,
    output:
        txt = resolve_results_filepath(
            "rseqc",
            "{sample}/se/{sample}.se.junction_annotation.txt"
        )
    params:
        bed = ref_path("rseqc", "bed"),
        out_prefix = resolve_results_filepath(
            "rseqc",
            "{sample}/se/{sample}.se"
        )
    benchmark:
        resolve_benchmarks_filepath(
            "rseqc",
            "{sample}.rseqc.junction_annotation.se.txt"
        )
    threads:
        conservative_cpu_count()
    conda:
        resolve_envs_filepath("rseqc.yaml")
    shell:
        "junction_annotation.py "
        "-r {params.bed} "
        "-i {input.bam} "
        "-o {params.out_prefix} "
        ">& {output.txt} "

rule rseqc_junction_saturation_pe:
    """
    RSeQC: junction saturation analysis for paired-end BAM.
    """
    input:
        bam = rules.star_align_pe.output.bam,
        bai = rules.index_bam_pe.output.bai,
    output:
        plot = resolve_results_filepath(
            "rseqc",
            "{sample}/{sample}.pe.junctionSaturation_plot.r"
        )
    params:
        bed = ref_path("rseqc", "bed"),
        out_prefix = resolve_results_filepath(
            "rseqc",
            "{sample}/{sample}.pe"
        )
    log:
        resolve_logs_filepath(
            "rseqc",
            "{sample}.rseqc.junction_saturation.pe.log"
        )
    benchmark:
        resolve_benchmarks_filepath(
            "rseqc",
            "{sample}.rseqc.junction_saturation.pe.txt"
        )
    threads:
        conservative_cpu_count()
    conda:
        resolve_envs_filepath("rseqc.yaml")
    shell:
        "junction_saturation.py "
        "-r {params.bed} "
        "-i {input.bam} "
        "-o {params.out_prefix} "
        ">& {log} "


rule rseqc_junction_saturation_se:
    """
    RSeQC: junction saturation analysis for single-end BAM.
    """
    input:
        bam = rules.star_align_se.output.bam,
        bai = rules.index_bam_se.output.bai,
    output:
        plot = resolve_results_filepath(
            "rseqc",
            "{sample}/{sample}.se.junctionSaturation_plot.r"
        )
    params:
        bed = ref_path("rseqc", "bed"),
        out_prefix = resolve_results_filepath(
            "rseqc",
            "{sample}/{sample}.se"
        )
    log:
        resolve_logs_filepath(
            "rseqc",
            "{sample}.rseqc.junction_saturation.se.log"
        )
    benchmark:
        resolve_benchmarks_filepath(
            "rseqc",
            "{sample}.rseqc.junction_saturation.se.txt"
        )
    threads:
        conservative_cpu_count()
    conda:
        resolve_envs_filepath("rseqc.yaml")
    shell:
        "junction_saturation.py "
        "-r {params.bed} "
        "-i {input.bam} "
        "-o {params.out_prefix} "
        ">& {log} "


rule rseqc_read_gc_pe:
    """
    RSeQC: GC content profile for paired-end BAM.
    """
    input:
        bam = rules.star_align_pe.output.bam,
        bai = rules.index_bam_pe.output.bai
    output:
        xls = resolve_results_filepath(
            "rseqc",
            "{sample}/{sample}.pe.GC.xls"
        )
    params:
        out_prefix = resolve_results_filepath(
            "rseqc",
            "{sample}/{sample}.pe"
        )
    benchmark:
        resolve_benchmarks_filepath(
            "rseqc",
            "{sample}.rseqc.read_GC.pe.txt"
        )
    threads:
        conservative_cpu_count()
    conda:
        resolve_envs_filepath("rseqc.yaml")
    shell:
        "read_GC.py "
        "-i {input.bam} "
        "-o {params.out_prefix} "


rule rseqc_read_gc_se:
    """
    RSeQC: GC content profile for seungle-end BAM.
    """
    input:
        bam = rules.star_align_se.output.bam,
        bai = rules.index_bam_se.output.bai
    output:
        xls = resolve_results_filepath(
            "rseqc",
            "{sample}/se/{sample}.se.GC.xls"
        )
    params:
        out_prefix = resolve_results_filepath(
            "rseqc",
            "{sample}/se/{sample}.se"
        )
    benchmark:
        resolve_benchmarks_filepath(
            "rseqc",
            "{sample}.rseqc.read_GC.se.txt"
        )
    threads:
        conservative_cpu_count()
    conda:
        resolve_envs_filepath("rseqc.yaml")
    shell:
        "read_GC.py "
        "-i {input.bam} "
        "-o {params.out_prefix} "

rule rseqc_read_duplication_pe:
    """
    RSeQC: read duplication profile for paired-end BAM.
    """
    input:
        bam = rules.star_align_pe.output.bam,
        bai = rules.index_bam_pe.output.bai
    output:
        xls = resolve_results_filepath(
            "rseqc",
            "{sample}/{sample}.pe.pos.DupRate.xls"
        )
    params:
        out_prefix = resolve_results_filepath(
            "rseqc",
            "{sample}/{sample}.pe"
        )
    benchmark:
        resolve_benchmarks_filepath(
            "rseqc",
            "{sample}.rseqc.read_duplication.pe.txt"
        )
    threads:
        conservative_cpu_count()
    conda:
        resolve_envs_filepath("rseqc.yaml")
    shell:
        "read_duplication.py "
        "-i {input.bam} "
        "-o {params.out_prefix} "

rule rseqc_read_duplication_se:
    """
    RSeQC: read duplication profile for sinle-end BAM.
    """
    input:
        bam = rules.star_align_se.output.bam,
        bai = rules.index_bam_se.output.bai
    output:
        xls = resolve_results_filepath(
            "rseqc",
            "{sample}/se/{sample}.se.pos.DupRate.xls"
        )
    params:
        out_prefix = resolve_results_filepath(
            "rseqc",
            "{sample}/se/{sample}.se"
        )
    benchmark:
        resolve_benchmarks_filepath(
            "rseqc",
            "{sample}.rseqc.read_duplication.se.txt"
        )
    threads:
        conservative_cpu_count()
    conda:
        resolve_envs_filepath("rseqc.yaml")
    shell:
        "read_duplication.py "
        "-i {input.bam} "
        "-o {params.out_prefix} "


rule rseqc_rpkm_saturation_pe:
    """
    RSeQC: RPKM saturation analysis for paired-end BAM.
    """
    input:
        bam = rules.star_align_pe.output.bam,
        bai = rules.index_bam_pe.output.bai,
    output:
        pdf = resolve_results_filepath(
            "rseqc",
            "{sample}/{sample}.pe.saturation.pdf"
        )
    params:
        bed = ref_path("rseqc", "bed"),
        out_prefix = resolve_results_filepath(
            "rseqc",
            "{sample}/{sample}.pe"
        )
    log:
        resolve_logs_filepath(
            "rseqc",
            "{sample}.rseqc.RPKM_saturation.pe.log"
        )
    benchmark:
        resolve_benchmarks_filepath(
            "rseqc",
            "{sample}.rseqc.RPKM_saturation.pe.txt"
        )
    threads:
        conservative_cpu_count()
    conda:
        resolve_envs_filepath("rseqc.yaml")
    shell:
        "RPKM_saturation.py "
        "-r {input.bed} "
        "-i {input.bam} "
        "-o {params.out_prefix} "
        ">& {log} "


rule rseqc_rpkm_saturation_se:
    """
    RSeQC: RPKM saturation analysis for single-end BAM.
    """
    input:
        bam = rules.star_align_se.output.bam,
        bai = rules.index_bam_se.output.bai,
    output:
        pdf = resolve_results_filepath(
            "rseqc",
            "{sample}/{sample}.se.saturation.pdf"
        )
    params:
        bed = ref_path("rseqc", "bed"),
        out_prefix = resolve_results_filepath(
            "rseqc",
            "{sample}/se/{sample}.se"
        )
    log:
        resolve_logs_filepath(
            "rseqc",
            "{sample}.rseqc.RPKM_saturation.se.log"
        )
    benchmark:
        resolve_benchmarks_filepath(
            "rseqc",
            "{sample}.rseqc.RPKM_saturation.se.txt"
        )
    threads:
        conservative_cpu_count()
    conda:
        resolve_envs_filepath("rseqc.yaml")
    shell:
        "RPKM_saturation.py "
        "-r {input.bed} "
        "-i {input.bam} "
        "-o {params.out_prefix} "
        ">& {log} "
