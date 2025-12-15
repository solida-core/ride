# RIDE – RNA-Seq Differential Expression Pipeline

RIDE is a reproducible and modular Snakemake pipeline for RNA-Seq differential expression analysis.  
It performs trimming, alignment, quantification, QC, and DE analysis using STAR, kallisto, and DESeq2.  
Each run is fully isolated and includes a snapshot of the configuration for complete traceability.

---

## 📦 Requirements

- Conda / Miniconda
- Snakemake ≥ 8.x
- Bash shell (Linux/macOS)

---

## 📥 Clone the Repository

Before installing the environment, clone the pipeline:

```bash
git clone https://github.com/solida-core/ride.git
cd ride
````

---

## 🔧 Setup

Install the Conda environment:

```bash
make install
```
By default, the environment is created in the standard Conda locations (eg: $HOME/anaconda/envs/ride)
You can specify a custom installation prefix using the prefix parameter:

```bash
make install prefix=/path/to/conda/envs/rna_seq_de_analysis
```

To update the environment:

```bash
make update
```

Show all available Make targets:

```bash
make help
```

---

## 🚀 Running the Pipeline

Launch the workflow with:

```bash
./run_ride.sh [options]
```

By default, each execution creates:

```
runs/<timestamp>/
```

### Available options

```
Usage: run_ride.sh [-h] [-n] [-s SNAKEFILE] [-c CONFIG_FILE] [-w WORKDIR] [-P PROFILE] [-p "PARAMS"]
Launch the RIDE Snakemake workflow.

Options:
-h                Show this help message and exit
-n                Run in dry-run mode (equivalent to --dry-run)
-s SNAKEFILE      Snakefile to use (default: workflow/Snakefile)
-c CONFIG_FILE    Configuration file (default: config/config.yaml)

-w WORKDIR        Run directory:
                    • empty  → runs/<timestamp>/
                    • name   → runs/<name>/
                    • path   → used exactly as provided

-P PROFILE        Snakemake profile directory (e.g. profiles/slurm)
-p PARAMS         Additional Snakemake parameters (passed verbatim)

```

### Examples

```bash
# Default run (creates runs/<timestamp>/)
./run_ride.sh

# Custom run name
./run_ride.sh -w test_run

# Use a different configuration file
./run_ride.sh -c config/hg38.yaml -w experiment_A

# Run with a Snakemake cluster profile
./run_ride.sh -P slurm_profile

# Pass extra Snakemake options
./run_ride.sh -p "--cores 40 --latency-wait 60"
````

---

## 🗂️ Project Structure

```
config/            # Configuration files
workflow/          # Snakefile, rules, scripts, envs
run_ride.sh        # Launcher
Makefile           # Environment automation
environment.yaml   # Conda environment
runs/              # Output run directories
```

Each run directory contains a copy of the configuration and all workflow outputs:

* logs
* QC reports
* STAR / kallisto outputs
* DESeq2 results
* `summary.tsv`
* `report.html`

---

## 📊 Output Summary

RIDE produces:

* Trimmed reads and QC reports
* STAR alignments
* kallisto quantification
* DESeq2 differential expression results
* Snakemake HTML report
* Execution summary table

---

## 🛠️ Makefile Targets

```makefile
make help        # Show help
make install     # Install the Conda environment
make update      # Recreate the environment
make clean       # Remove the environment
make run         # Execute the pipeline via the launcher
```

---

## 👤 Author

Maintained by Rossano Atzeni — CRS4 Bioinformatics Unit

---

## 📄 License

Distributed under the MIT License (see `LICENSE`).
