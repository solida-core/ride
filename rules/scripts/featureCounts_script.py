from subprocess import run

if len(snakemake.input[0]) > 1:
    cmd = [snakemake.params['cmd'],
           '-p',
           '-T', str(snakemake.threads),
           '-t', str(snakemake.params['gtf_feature_type']),
           '-a', str(snakemake.params['gtf']),
           '-o', snakemake.output,
           snakemake.input[1]]
else:
    cmd = [snakemake.params['cmd'],
           '-T', str(snakemake.threads),
           '-t', str(snakemake.params['gtf_feature_type']),
           '-a', str(snakemake.params['gtf']),
           '-o', snakemake.output,
           snakemake.input[1]]
run(cmd)
