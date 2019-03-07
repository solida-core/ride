rule rseqc_bam_stat:
    input:
        "star/{sample}/{sample}.Aligned.sortedByCoord.out.bam"
    output:
        "rseqc/{sample}/{sample}.bam_stat.txt"
    conda:
        "../envs/rseqc.yaml"
    shell:
        "bam_stat.py "
        "-i {input} "
        ">& {output}"

rule rseqc_genebody_coverage:
    input:
        bam="star/{sample}/{sample}.Aligned.sortedByCoord.out.bam",
        bai="star/{sample}/{sample}.Aligned.sortedByCoord.out.bam.bai"
    output:
        "rseqc/{sample}/{sample}.geneBodyCoverage.txt"
    params:
        out_basename="rseqc/{sample}/{sample}",
#        ref="/ELS/els9/users/biosciences/references/rseqc/hg19.HouseKeepingGenes.bed"
        ref=resolve_single_filepath(*references_abs_path(ref="rseqc_reference"), config.get("housekeeping"))
    conda:
        "../envs/rseqc.yaml"
    log:
        "logs/rseqc/genebody_coverage/{sample}_genebodycoverage.log"
    shell:
        "geneBody_coverage.py "
        "-r {params.ref} "
        "-i {input.bam} "
        "-o {params.out_basename} "
        ">& {log} "

rule rseqc_junction_annotation:
    input:
        bam="star/{sample}/{sample}.Aligned.sortedByCoord.out.bam",
        bai="star/{sample}/{sample}.Aligned.sortedByCoord.out.bam.bai"
    output:
        out="rseqc/{sample}/{sample}.junction.txt"
    params:
        out_basename="rseqc/{sample}/{sample}",
        ref=resolve_single_filepath(*references_abs_path(ref="rseqc_reference"), config.get("refseq"))
    conda:
        "../envs/rseqc.yaml"
    shell:
        "junction_annotation.py "
        "-r {params.ref} "
        "-i {input.bam} "
        "-o {params.out_basename} "
        ">& {output.out}"

rule rseqc_junction_saturation:
    input:
        bam="star/{sample}/{sample}.Aligned.sortedByCoord.out.bam",
        bai="star/{sample}/{sample}.Aligned.sortedByCoord.out.bam.bai"
    output:
        plotr="rseqc/{sample}/{sample}.junctionSaturation_plot.r"
    params:
        out_basename="rseqc/{sample}/{sample}",
        ref=resolve_single_filepath(*references_abs_path(ref="rseqc_reference"), config.get("refseq"))
    conda:
        "../envs/rseqc.yaml"
    log:
        "rseqc/{sample}/{sample}.junctionSaturation.txt"
    shell:
        "junction_saturation.py "
        "-r {params.ref} "
        "-i {input.bam} "
        "-o {params.out_basename} "
        "&> {log}"


rule rseqc_GC:
    input:
        bam="star/{sample}/{sample}.Aligned.sortedByCoord.out.bam",
        bai="star/{sample}/{sample}.Aligned.sortedByCoord.out.bam.bai"
    output:
        "rseqc/{sample}/{sample}.GC.xls"
    params:
        out_basename="rseqc/{sample}/{sample}"
    conda:
        "../envs/rseqc.yaml"
    shell:
        "read_GC.py "
        "-i {input.bam} "
        "-o {params.out_basename}"

rule rseqc_read_distribution:
    input:
        bam="star/{sample}/{sample}.Aligned.sortedByCoord.out.bam",
        bai="star/{sample}/{sample}.Aligned.sortedByCoord.out.bam.bai"
    output:
        "rseqc/{sample}/{sample}.read_distribution.txt"
    params:
        ref=resolve_single_filepath(*references_abs_path(ref="rseqc_reference"), config.get("refseq"))
    conda:
        "../envs/rseqc.yaml"
    shell:
        "read_distribution.py "
        "-r {params.ref} "
        "-i {input.bam} "
        ">& {output}"


rule rseqc_infer_experiment:
    input:
        bam="star/{sample}/{sample}.Aligned.sortedByCoord.out.bam",
        bai="star/{sample}/{sample}.Aligned.sortedByCoord.out.bam.bai"
    output:
        "rseqc/{sample}/{sample}.infer_experiment.txt"
    params:
        ref=resolve_single_filepath(*references_abs_path(ref="rseqc_reference"), config.get("refseq"))
    conda:
        "../envs/rseqc.yaml"
    shell:
        "infer_experiment.py "
        "-r {params.ref} "
        "-i {input.bam} "
        ">& {output}"


rule rseqc_read_duplication:
    input:
        bam="star/{sample}/{sample}.Aligned.sortedByCoord.out.bam",
        bai="star/{sample}/{sample}.Aligned.sortedByCoord.out.bam.bai"
    output:
        out1="rseqc/{sample}/{sample}.pos.DupRate.xls"
    params:
        out_basename="rseqc/{sample}/{sample}"
    conda:
        "../envs/rseqc.yaml"
    shell:
        "read_duplication.py "
        "-i {input.bam} "
        "-o {params.out_basename}"

rule rseqc_RPKM_saturation:
    input:
        bam="star/{sample}/{sample}.Aligned.sortedByCoord.out.bam",
        bai="star/{sample}/{sample}.Aligned.sortedByCoord.out.bam.bai"
    output:
        out1="rseqc/{sample}/{sample}.saturation.pdf"
    params:
        out_basename="rseqc/{sample}/{sample}",
        ref=resolve_single_filepath(*references_abs_path(ref="rseqc_reference"), config.get("refseq"))
    log:
        "logs/rseqc/{sample}.RPKM_saturation.log"
    conda:
        "../envs/rseqc.yaml"
    shell:
        "RPKM_saturation.py "
        "-r {params.ref} "
        "-i {input.bam} "
        "-o {params.out_basename} "
        ">& {log} "

