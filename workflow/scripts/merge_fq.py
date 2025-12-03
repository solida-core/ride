from subprocess import run
import shutil
import os

inputs = list(snakemake.input)
output = snakemake.output.fq
outdir = snakemake.output.dir
logfile = snakemake.log[0] if snakemake.log else None

# Ensure output directory exists
os.makedirs(outdir, exist_ok=True)

# Ensure log directory exists
if logfile:
    os.makedirs(os.path.dirname(logfile), exist_ok=True)

def log(msg):
    """Write msg to log file."""
    if logfile:
        with open(logfile, "a") as lf:
            lf.write(msg + "\n")

log("=== FASTQ MERGE START ===")
log(f"Output: {output}")
log(f"Input files ({len(inputs)}):")
for f in inputs:
    log(f" - {f}")

if len(inputs) > 1:
    cmd = ["cat"] + inputs
    log("Merging files via: " + " ".join(cmd))
    with open(output, "wb") as fout:
        run(cmd, stdout=fout, check=True)
    log("Merge completed.")
else:
    log("Single input file → copying.")
    shutil.copyfile(inputs[0], output)
    log("Copy completed.")

log("=== FASTQ MERGE END ===")

