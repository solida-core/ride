from subprocess import run
import shutil

inputs = list(snakemake.input)
output = snakemake.output[0]

print(inputs)
print(len(inputs))

if len(inputs) > 1:
    # Multiple input files → concatenate them using `cat`
    # Use argument list instead of shell=True for safety
    cmd = ["cat"] + inputs
    with open(output, "wb") as fout:
        run(cmd, stdout=fout, check=True)
else:
    # Single input file → just copy it
    shutil.copyfile(inputs[0], output)
