# ============================================================
# RiDE - RNA-Seq Differential Expression Snakemake Pipeline
# Makefile
# ============================================================

SHELL=/bin/bash

ENV_NAME=ride
ENV_FILE=environment.yaml
ENV_METADATA=.ride_envpath
SCRIPT_PATH=run_ride.sh

# User can override install path:  make install prefix=/path/to/env
prefix ?=

.PHONY: help check install update clean run

# ------------------------------------------------------------
# Help
# ------------------------------------------------------------
help:
	@echo ""
	@echo "Usage:"
	@echo "  make install [prefix=/path]   Create Conda environment"
	@echo "  make update [prefix=/path]    Recreate environment"
	@echo "  make clean                    Remove environment"
	@echo "  make run                      Execute the RIDE pipeline"
	@echo ""
	@echo "Run parameters:"
	@echo "    w=<workdir>   Run directory"
	@echo "    c=<config>    Config file"
	@echo "    P=<profile>   Snakemake profile"
	@echo "    p=\"params\"   Extra Snakemake parameters"
	@echo ""
	@echo "Examples:"
	@echo "    make install"
	@echo "    make install prefix=/opt/envs/ride"
	@echo "    make run w=test_run"
	@echo ""

# ------------------------------------------------------------
# Check for Conda installation
# ------------------------------------------------------------
check:
	@echo "Checking for Conda..."
	@which conda >/dev/null 2>&1 || { \
		echo "ERROR: Conda not found. Install Miniconda or Anaconda."; \
		exit 1; \
	}
	@echo "Conda detected."

# ------------------------------------------------------------
# Install environment (with optional --prefix)
# ------------------------------------------------------------
install: check
	@echo "Installing environment '$(ENV_NAME)'"

	@if [ -n "$(prefix)" ]; then \
		ENV_PREFIX="$(prefix)"; \
	else \
		ENV_PREFIX="$$HOME/.conda/envs/$(ENV_NAME)"; \
	fi; \
	echo "Environment prefix: $$ENV_PREFIX"; \
	\
	if [ -d "$$ENV_PREFIX" ]; then \
		echo "Environment already exists at $$ENV_PREFIX"; \
		read -p "Overwrite it? [Y/n] " ans; \
		if ! echo $$ans | grep -iq "^y" ; then \
			echo "Skipping installation."; exit 0; \
		fi; \
		echo "Removing old environment..."; \
		conda remove --prefix $$ENV_PREFIX --all -y; \
	fi; \
	\
	echo "Creating new environment..."; \
	conda env create --prefix $$ENV_PREFIX -f $(ENV_FILE); \
	echo $$ENV_PREFIX > $(ENV_METADATA); \
	echo "Environment path saved to $(ENV_METADATA)"; \
	chmod +x $(SCRIPT_PATH); \
	echo "Installation complete."

# ------------------------------------------------------------
# Update environment (always recreates)
# ------------------------------------------------------------
update: check
	@if [ -n "$(prefix)" ]; then \
		ENV_PREFIX="$(prefix)"; \
	elif [ -f $(ENV_METADATA) ]; then \
		ENV_PREFIX=$$(cat $(ENV_METADATA)); \
	else \
		ENV_PREFIX="$$HOME/.conda/envs/$(ENV_NAME)"; \
	fi; \
	echo "Recreating environment at $$ENV_PREFIX"; \
	conda remove --prefix $$ENV_PREFIX --all -y || true; \
	conda env create --prefix $$ENV_PREFIX -f $(ENV_FILE); \
	echo $$ENV_PREFIX > $(ENV_METADATA); \
	chmod +x $(SCRIPT_PATH); \
	echo "Environment recreated."

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
	@if [ -f $(ENV_METADATA) ]; then \
		ENV_PREFIX=$$(cat $(ENV_METADATA)); \
		echo "Removing environment at $$ENV_PREFIX"; \
		conda remove --prefix $$ENV_PREFIX --all -y || true; \
		rm -f $(ENV_METADATA); \
	else \
		echo "No environment metadata found. Nothing to clean."; \
	fi
	@echo "Environment removed."

