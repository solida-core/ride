rule star_build_index:
    input:
        resolve_single_filepath(*references_abs_path(ref='genome_reference'),
                                config.get("genome_fasta"))
    output:
          genomeDir='star/index',
          length='star/index/chrLength.txt'
    conda:
        "envs/star2.yaml"
    params:
        gtf=resolve_single_filepath(*references_abs_path(ref='genes_reference'),
                                    config.get("genes_gtf"))
    threads: pipeline_cpu_count()
    shell:
          "STAR "
          "--runMode genomeGenerate "
          "--runThreadN {threads} "
          "--genomeDir {output.genomeDir} "
          "--genomeFastaFiles {input} "
          "--sjdbGTFfile {params.gtf} "
          "--sjdbOverhang 100"


rule star_map:
    input:
        lambda wildcards: config["samples"][wildcards.sample],
        length=rules.star_build_index.output.length
    output:
        out1="star/{sample}/{sample}.Aligned.sortedByCoord.out.bam",

    conda:
        "envs/star2.yaml"
    params:
        genomedir = 'star/index',
        sample = "de",
        platform='platform',
        center='center',
        out_basename="star/{sample}/{sample}."
    threads: pipeline_cpu_count()
    log:
        "logs/star/{sample}/{sample}_star_map.log"
    shell:
        "STAR "
        "--runMode alignReads "
        "--genomeDir {params.genomedir} "
        r" --outSAMattrRGline  ID:{params.sample} SM:{params.sample} PL:{params.platform}  PU:{params.platform} CN:{params.center} "
        "--readFilesIn {input} "
        "--readFilesCommand zcat "
        "--outStd Log "
        "--outSAMunmapped Within "
        "--outSAMtype BAM SortedByCoordinate "
        "--outWigType wiggle "
        "--outWigStrand Stranded "
        "--runThreadN {threads} "
        "--outFileNamePrefix {params.out_basename} "
        "2> {log} "

rule featureCounts_run:
    input:
        lambda wildcards: config["samples"][wildcards.sample],
        bam="star/{sample}/{sample}.Aligned.sortedByCoord.out.bam"
    output:
        "star/{sample}/count/{sample}_counts.cnt"
    conda:
        "envs/featureCounts.yaml"
    params:
         cmd="featureCounts",
         gtf=resolve_single_filepath(*references_abs_path(ref='genes_reference'),
                                    config.get("genes_gtf")),
         gtf_feature_type=config.get("rules").get("featureCounts_run").get("gtf_feature_type"),
         gtf_attribute_type=config.get("rules").get("featureCounts_run").get("gtf_attribute_type"),
    threads: pipeline_cpu_count()
    script:
         "scripts/featureCounts_script.py"


rule samtools_index:
    input:
        "star/{sample}/{sample}.Aligned.sortedByCoord.out.bam"
    output:
        "star/{sample}/{sample}.Aligned.sortedByCoord.out.bam.bai"
    conda:
        "envs/samtools.yaml"
    benchmark:
        "benchmarks/samtools/index/{sample}.txt"
    shell:
        "samtools index "
        "{input} "
