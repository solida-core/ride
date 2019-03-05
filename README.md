# ride - RNA seq data analysis pipeline
Snakemake pipeline for RNA-seq data analysis.

## Workflow
ride workflow comprise two analysis phases:
 * _concatenate_, from multiple-fastq to single fastq  
 * _ride_, from single fastq to final analysis

## PRE-ANALYSIS STEP
In a pre-analysis step, multiple-fastq files have to be merged in a single fastq for sample. To do so is provided a script "concatenate.py" which manages this step. The script also create a config.file specific for the project with replaced samples and paths.
The script requires:
 * input file: a file yaml including the complete paths for fastqs for sample
 * folder: a folder where merged fastq files will be stored. Note that this folder mustn't exist (to avoid overwriting), otherwise the script will raise an error.
 * project_name: the project name. It will be included in the generated config filename (default is project folder name).
```bash
python concatenate.py -h
python concatenate.py --input_file INPUT_FILE --folder PATH [--project_name PATH]
``` 
## Analysis
```bash
./run.project.sh -s Snakefile -c config.[project_name].yaml
```

