
rule star_build_index:
    input:
        resolve_single_filepath(*references_abs_path(), config.get("genome_fasta"))
    output:
        length='star/index/chrLength.txt'
    conda:
        "../envs/star2.yaml"
    params:
        gtf=resolve_single_filepath(*references_abs_path(ref='references'), config.get("genes_gtf")),
        genomeDir='star/index'
    threads: pipeline_cpu_count()
    shell:
        "STAR "
        "--runMode genomeGenerate "
        "--runThreadN {threads} "
        "--genomeDir {params.genomeDir} "
        "--genomeFastaFiles {input} "
        "--sjdbGTFfile {params.gtf} "
        "--sjdbOverhang 100"



rule star_map:
    input:
        "reads/trimmed/{sample}-R1-trimmed.fq.gz",
        "reads/trimmed/{sample}-R2-trimmed.fq.gz",
        length=rules.star_build_index.output.length
    output:
        out1="star/{sample}/{sample}.Aligned.sortedByCoord.out.bam",
        out2="star/{sample}/{sample}.Log.final.out"
    conda:
        "../envs/star2.yaml"
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
        "--readFilesIn {input[0]} {input[1]} "
        "--readFilesCommand zcat "
        "--outStd Log "
        "--outSAMunmapped Within "
        "--outSAMtype BAM SortedByCoordinate "
        "--outWigType wiggle "
        "--outWigStrand Stranded "
        "--runThreadN {threads} "
        "--outFileNamePrefix {params.out_basename} "
        "2> {log} "



rule samtools_index:
    input:
        "star/{sample}/{sample}.Aligned.sortedByCoord.out.bam"
    output:
        "star/{sample}/{sample}.Aligned.sortedByCoord.out.bam.bai"
    conda:
        "../envs/samtools.yaml"
    benchmark:
        "benchmarks/samtools/index/{sample}.txt"
    shell:
        "samtools index "
        "{input} "
