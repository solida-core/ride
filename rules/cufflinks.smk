rule assembly:
    input:
        expand('mapped/{sample}/accepted_hits.bam', sample=config.get("samples"))
    output:
        'assembly/{sample}/transcripts.gtf',
        dir='assembly/{sample}'
    conda:
        "../envs/cufflinks.yaml"
    params:
        genome="{label}.fa".format(label=get_references_label())
    threads: pipeline_cpu_count()
    shell:
        'cufflinks --num-threads {threads} -o {output.dir} '
        '--frag-bias-correct {params.genome} {input}'

rule compose_merge:
    input:
        expand('assembly/{sample}/transcripts.gtf', sample=config.get("samples"))
    output:
        txt='assembly/assemblies.txt'
    run:
        with open(output.txt, 'w') as out:
            print(*input, sep="\n", file=out)

rule merge_assemblies:
    input:
        'assembly/assemblies.txt'
    output:
        'assembly/merged/merged.gtf', dir='assembly/merged'
    conda:
        "../envs/cufflinks.yaml"
    params:
        genome="{label}.fa".format(label=get_references_label())
    shell:
        'cuffmerge -o {output.dir} -s {params.genome} {input}'

rule compare_assemblies:
    input:
        'assembly/merged/merged.gtf'
    output:
        'assembly/comparison/all.stats',
        dir='assembly/comparison'
    conda:
        "../envs/cufflinks.yaml"
    params:
        genome="{label}.fa".format(label=get_references_label()),
        gtf=resolve_single_filepath(*references_abs_path(ref='genes_reference'),
                                    config.get("genes_gtf"))
    shell:
        'cuffcompare -o {output.dir}/all -s {params.genome} -r {params.gtf} '
        '{input}'

rule diffexp:
    input:
        class1=expand('mapped/{sample}/accepted_hits.bam', sample=config.get('classes').get('C1')),
        class2=expand('mapped/{sample}/accepted_hits.bam', sample=config.get('classes').get('C2')),
        gtf='assembly/merged/merged.gtf'
    output:
        'diffexp/gene_exp.diff', 'diffexp/isoform_exp.diff'
    params:
        class1=lambda wildcards, input: ",".join(input.class1),
        class2=lambda wildcards, input: ",".join(input.class2)
    threads: pipeline_cpu_count()
    shell:
        "cuffdiff "
        "--num-threads {threads} "
        "{input.gtf} "
        "{params.class1} "
        "{params.class2}"