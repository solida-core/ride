
## kallisto abundance tsv nella cartella abundance/sample.tsv
## star_bam nella cartella bams/sample.bam sample.bam.bai
## qc multiqc.html

def get_sample_by_client(wildcards, reheader, label='LIMS', structure="folder/{sample}.extension"):
    re.sub(r"{sample}",reheader.loc[wildcards.Client,[label]][0], structure)
    return re.sub(r"{sample}",reheader.loc[wildcards.Client,[label]][0], structure)

rule delivery_kallisto:
    input:
        lambda wildcards: get_sample_by_client(wildcards, reheader, label=config.get("internal_sid"), structure='kallisto/{sample}/abundance.tsv')
    output:
        "delivery/abundance/{Client}.tsv"
    shell:
        "cp {input} {output}"


rule delivery_star:
    input:
        bam=lambda wildcards: get_sample_by_client(wildcards, reheader, label=config.get("internal_sid"), structure='star/{sample}/{sample}.Aligned.sortedByCoord.out.bam'),
        bai=lambda wildcards: get_sample_by_client(wildcards, reheader, label=config.get("internal_sid"), structure='star/{sample}/{sample}.Aligned.sortedByCoord.out.bam.bai')
    output:
        bam="delivery/bams/{Client}.bam",
        bai="delivery/bams/{Client}.bam.bai"
    shell:
        "cp {input.bam} {output.bam} && "
        "cp {input.bai} {output.bai} "


rule delivery_multiqc:
    input:
        "qc/multiqc.html"
    output:
        "delivery/qc/multiqc.html"
    shell:
        "cp {input} {output}"