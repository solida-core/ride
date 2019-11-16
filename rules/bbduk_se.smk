rule bbduk_se:
    input:
       reads="reads/trimmed/{sample}-R1-trimmed.fq.gz"
    output:
       nonribo="reads/bbduk/{sample}.fa",
       stats="logs/bbduk/{sample}.log"
    params:
       ref=resolve_single_filepath(*references_abs_path(ref="bbduk_reference"), config.get("bbduk_se"))
    shell:
       "bbduk.sh "
       "in={input} "
       "out={output.nonribo} "
       "ref={params.ref} "
       "stats={output.stats}"