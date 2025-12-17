# RIDE Configuration Files

This directory contains all configuration files required to run the RIDE pipeline.
Each run uses these files to define samples, input data, reference resources, and analysis parameters.

The full `config/` directory is copied into each run folder to guarantee complete reproducibility.

---

## 1. Main configuration file: `config.yaml`

The main configuration file defines input tables, reference resources, and tool-specific parameters.

### Input tables

```yaml
samples: config/samples.tsv
units: config/units.tsv
reheader: config/reheader.tsv
````

* `samples.tsv`: defines biological samples and their associated sequencing units
* `units.tsv`: defines FASTQ files for each sequencing unit
* `reheader.tsv`: optional mapping between internal sample IDs and human-readable names

---

### Reference resources

#### Genome reference

```yaml
resources:
  genome:
    basepath: "/path/to/refdata"
    provider: "ensembl"
    release: "GRCh38"
    fasta: "Homo_sapiens.GRCh38.dna.primary_assembly.fa"
    gtf: "Homo_sapiens.GRCh38.115.gtf"
```

* `basepath`: directory containing reference files
* `fasta`: genome FASTA file
* `gtf`: gene annotation file (GTF format)

#### Transcriptome reference (kallisto)

```yaml
transcriptome:
  basepath: "/path/to/refdata"
  provider: "ensembl"
  release: "GRCh38"
  cdna: "Homo_sapiens.GRCh38.cdna.all.fa.gz"
```

Used to build the kallisto index.

#### RSeQC resources

```yaml
rseqc:
  basepath: "/path/to/refdata"
  provider: "rseqc"
  release: "hg38"
  housekeeping: "hg38.HouseKeepingGenes.bed"
  refseq: "hg38_RefSeq.bed"
```

Used for RNA-seq quality control metrics.

---

### Trimming parameters

```yaml
trimming:
  quality: 20
  min_length: 20
  trim_poly_g: true
  trim_poly_x: true
  detect_adapter: true
```

Defines read trimming and adapter detection settings.

---

### Kallisto parameters

```yaml
kallisto:
  index_name: "kallisto_index.idx"
  frag_len: 280
  frag_sd: 80
  boot: 100
```

* `frag_len`, `frag_sd`: fragment length parameters (single-end data)
* `boot`: number of bootstrap samples

---

### STAR parameters

```yaml
star:
  index_dir: "index"
  sjdb_overhang: 100
```

* `index_dir`: STAR genome index directory
* `sjdb_overhang`: read length minus one

---

## 2. Sample definition: `samples.tsv`

Defines biological samples and their associated sequencing units.

```tsv
sample	units
ERS179576	HSQ1008_141.L005.ERS179576,HSQ1008_141.L007.ERS179576,HSQ1009_88.L001.ERS179576
ERS179577	HSQ1009_86.L001.ERS179577
```

* `sample`: unique biological sample identifier
* `units`: comma-separated list of sequencing units belonging to the sample

---

## 3. Sequencing units: `units.tsv`

Defines FASTQ files associated with each sequencing unit.

```tsv
sample	unit	fq1	fq2
ERS179576	HSQ1008_141.L005.ERS179576	path_to_datasets/ERR174310_1.fastq.gz	path_to_datasets/ERR174310_2.fastq.gz
ERS179576	HSQ1008_141.L007.ERS179576	path_to_datasets/ERR174312_1.fastq.gz	path_to_datasets/ERR174312_2.fastq.gz
ERS179576	HSQ1009_88.L001.ERS179576	path_to_datasets/ERR174314_1.fastq.gz
ERS179577	HSQ1009_86.L001.ERS179577	path_to_datasets/ERR174324_1.fastq.gz	path_to_datasets/ERR174324_2.fastq.gz
```

* `sample`: biological sample ID (must match `samples.tsv`)
* `unit`: sequencing unit / lane / library identifier
* `fq1`: FASTQ file (read 1)
* `fq2`: FASTQ file (read 2, optional for single-end data)

---

## 4. Sample renaming: `reheader.tsv`

Optional mapping from internal sample IDs to human-readable names used in reports and plots.

```tsv
sample_id	sample_name
ERS179576	76
ERS179577	77
```

* `sample_id`: original sample identifier
* `sample_name`: label used in downstream outputs

---

## Notes

* All paths can be absolute or relative
* Samples with missing FASTQ files will cause the pipeline to fail
* The configuration directory is copied into each run directory to ensure full traceability

````


