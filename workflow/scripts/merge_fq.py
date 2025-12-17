import os
import shutil
import sys
from subprocess import run

inputs = list(snakemake.input)
output = snakemake.output.fq
outdir = os.path.dirname(snakemake.output.fq)
logfile = snakemake.log[0] if snakemake.log else None

# === EXIT IMMEDIATELY IF NO INPUT FILES ===
if len(inputs) == 0:
    sys.exit(0)  # Exit cleanly without creating logs or outputs

# Ensure output directory exists
os.makedirs(outdir, exist_ok=True)

# Ensure log directory exists
if logfile:
    os.makedirs(os.path.dirname(logfile), exist_ok=True)


def log(msg):
    """Write msg to log file."""
    print(msg)
    if logfile:
        with open(logfile, "a") as lf:
            lf.write(msg + "\n")


log("=== FASTQ MERGE START ===")
log(f"Output: {output}")
log(f"Input files ({len(inputs)}):")
for f in inputs:
    log(f" - {f}")

if len(inputs) > 1:
    cmd = f"cat {' '.join(inputs)} > {output}"
    log("Merging files via: " + cmd)
    run(cmd, shell=True, check=True)
    log("Merge completed.")
else:
    log("Single input file → copying.")
    shutil.copyfile(inputs[0], output)
    log("Copy completed.")

log("=== FASTQ MERGE END ===")
