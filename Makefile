# ============================================================
# RIDE - RNA-Seq Differential Expression Snakemake Pipeline
# Makefile
# ============================================================

SHELL=/bin/bash

ENV_NAME=ride
ENV_FILE=environment.yaml
ENV_METADATA=.ride_envpath
SCRIPT_PATH=run_ride.sh

.PHONY: help check install update clean run detect-path

# ------------------------------------------------------------
# Help
# ------------------------------------------------------------
help:
	@echo ""
	@echo "Usage:"
	@echo "  make install        Create (or optionally overwrite) Conda environment '$(ENV_NAME)'"
	@echo "  make update         Force recreate the environment from $(ENV_FILE)"
	@echo "  make clean          Remove the Conda environment"
	@echo "  make run            Execute the RIDE pipeline (supports parameters)"
	@echo ""
	@echo "Run parameters (optional):"
	@echo "    w=<workdir>       Run output directory"
	@echo "    c=<config>        Config file (default: config/config.yaml)"
	@echo "    P=<profile>       Snakemake profile (e.g. profiles/slurm)"
	@echo "    p=\"<params>\"     Additional Snakemake parameters"
	@echo ""
	@echo "Examples:"
	@echo "    make run"
	@echo "    make run w=test_run"
	@echo "    make run w=run1 c=config/tumor.yaml"
	@echo "    make run w=run1 P=profiles/slurm"
	@echo "    make run w=run1 p=\"--cores 40 --latency-wait 90\""
	@echo ""

# ------------------------------------------------------------
# Check for conda and mamba installation
# ------------------------------------------------------------
check:
	@echo "Checking for Conda..."
	@which conda >/dev/null 2>&1 || { \
		echo "ERROR: Conda not found. Install Miniconda or Anaconda."; \
		exit 1; \
	}
	@echo "Conda detected."

	@echo "Checking for Mamba..."
	@if ! conda list -n base | grep -q mamba; then \
		echo "Installing mamba in base environment..."; \
		conda install -n base -c conda-forge mamba -y; \
	else \
		echo "Mamba already installed."; \
	fi

# ------------------------------------------------------------
# Install environment (interactive overwrite)
# ------------------------------------------------------------
install: check
	@echo "Checking if environment '$(ENV_NAME)' exists..."
	@if conda env list | awk '{print $$1}' | grep -qx "$(ENV_NAME)"; then \
		echo "Environment '$(ENV_NAME)' already exists."; \
		read -p "Overwrite it? [Y/n] " ans; \
		if [[ "$$ans" =~ ^[Yy]$$ ]]; then \
			echo "Removing existing environment..."; \
			conda remove -n $(ENV_NAME) --all -y; \
			echo "Creating environment from $(ENV_FILE)..."; \
			mamba env create -f $(ENV_FILE); \
		else \
			echo "Skipping environment creation."; \
		fi; \
	else \
		echo "Creating environment from $(ENV_FILE)..."; \
		mamba env create -f $(ENV_FILE); \
	fi

	@$(MAKE) detect-path
	@chmod +x $(SCRIPT_PATH)
	@echo "Installation complete."

# ------------------------------------------------------------
# Detect environment path and store it in .ride_envpath
# ------------------------------------------------------------
detect-path:
	@echo "Detecting actual environment path..."
	@conda env list | awk -v env="$(ENV_NAME)" '\
		$$(NF)==env || $$1==env {print $$NF}' > $(ENV_METADATA)

	@if [[ ! -s $(ENV_METADATA) ]]; then \
		echo "WARNING: Could not detect environment path for $(ENV_NAME)."; \
		echo "Conda may not list it in a standard format."; \
	else \
		echo "Environment path stored in $(ENV_METADATA):"; \
		cat $(ENV_METADATA); \
	fi

# ------------------------------------------------------------
# Force recreation of the environment
# ------------------------------------------------------------
update: check
	@echo "Recreating environment '$(ENV_NAME)'..."
	@conda remove -n $(ENV_NAME) --all -y || true
	@mamba env create -f $(ENV_FILE)
	@$(MAKE) detect-path
	@chmod +x $(SCRIPT_PATH)
	@echo "Environment recreated."

# ------------------------------------------------------------
# Run the RIDE pipeline
# ------------------------------------------------------------
run:
	@echo "Executing RIDE pipeline"
	@./$(SCRIPT_PATH) \
	    ${w:+-w $(w)} \
	    ${c:+-c $(c)} \
	    ${P:+-P $(P)} \
	    ${p:+-p "$(p)"}

# ------------------------------------------------------------
# Clean/remove Conda environment
# ------------------------------------------------------------
clean:
	@echo "Removing conda environment '$(ENV_NAME)'..."
	@conda env remove --name $(ENV_NAME) -y || true
	@rm -f $(ENV_METADATA)
	@echo "Environment removed."
