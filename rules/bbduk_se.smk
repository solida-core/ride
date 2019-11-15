rule bbduk_se:
    input:
       reads="reads/trimmed/{sample}-R1-trimmed.fq.gz"
    output:
       nonribo="reads/bbduk/{sample}.fa",
       stats="logs/bbduk/{sample}.log"
    params:
       ref="config.get("rules").get("bbduk_se").get("ref")"
    shell:
       "bbduk.sh "
       "in={input} "
       "out={output.nonribo} "
       "ref={ref} "
       "stats={output.stats}"