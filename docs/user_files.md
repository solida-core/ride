## User files
**RiDE** (RNA Differential Expression) is a pipeline for **RNA-seq** data analysis.
The standardization provided by [solida-core]() requires perhaps accessory user-defined files to have a given organization.

#### Required Input Files
**RiDE** pipeline requires a single mandatory user input:
* [samples](#samplestsv)

This file contains **tab-separated** information about samples to be analyzed and must be declared into the `config.yaml` file:
```
samples: "path_to_input_files/samples.tsv"
```

#### samples.tsv
Samples file include information about sample IDs the absolute path of sequencing data file (or files if paired-end). 

The file structure is described in the example below::
```
sample  	fq1 	                                fq2
ERS179576	path_to_datasets/ERR174310_1.fastq.gz	path_to_datasets/ERR174310_2.fastq.gz
ERS179577	path_to_datasets/ERR174324_1.fastq.gz	path_to_datasets/ERR174324_2.fastq.gz
``` 
Where:

* sample: sample IDs;
* fq1-fq2: complete path and R1- and R2-fastq files (for paired-end sequencing).
______________________________________
