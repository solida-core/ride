
rule all:
    input:
        # Tophat+Cufflinks workflow
        'assembly/comparison/all.stats',
        'diffexp/isoform_exp.diff',
        # Kallisto+Sleuth workflow
        "kallisto/DEGS/gene_table.txt",
        # Star+htseq workflow
        expand("star/{sample}/count/{sample}_counts.cnt", sample=config.get('samples')),



include_prefix="rules"

include:
    include_prefix + "/notify.rules"
include:
    include_prefix + "/functions.py"
include:
    include_prefix + "/kallisto.rules"
include:
    include_prefix + "/tophat.rules"
include:
    include_prefix + "/cufflinks.rules"
include:
    include_prefix + "/star2.rules"
