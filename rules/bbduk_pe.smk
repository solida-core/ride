rule bbduk_pe:
     input:
         r1="reads/trimmed/{sample}-R1-trimmed.fq.gz",
         r2="reads/trimmed/{sample}-R2-trimmed.fq.gz"
     output:
         r1="reads/bbduk/r1_{sample}.fa"
         r2="reads/bbduk/r2_{sample}.fa"
         stats="logs/bbduk/{sample}.log"
     params:
        ref=resolve_single_filepath(*references_abs_path(ref="bbduk_reference"), config.get("bbduk_pe"))
     shell:
        "bbduk.sh "
        "in={input.r1} "
        "in2={input.r2} "
        "out={output.r1} "
        "out2={output.r2} "
        "ref={params.ref} "
        "stats={output.stats}"