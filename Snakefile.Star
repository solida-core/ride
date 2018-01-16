
rule all:
    input:
        # Star+htseq workflow
        expand("star/{sample}/count/{sample}_counts.cnt", sample=config.get
        ('samples')),



include_prefix="rules"

include:
    include_prefix + "/notify.rules"
include:
    include_prefix + "/functions.py"
include:
    include_prefix + "/star2.rules"
