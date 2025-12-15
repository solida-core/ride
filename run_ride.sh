#!/usr/bin/env bash

# ====================================================
# RiDE Snakemake launcher
#
# author: "Rossano Atzeni"
# pipeline: "RIDE RNA-Seq"
# version: "1.0.0"
# description: "Differential expression workflow"
# ====================================================

set -euo pipefail

# ========================
# User-editable settings
# ========================

# Directory where all runs will be stored (instead of hardcoded "results")
RUNS_DIR="runs"

# Default config file
DEFAULT_CONFIG_FILE="config/config.yaml"

SNAKE_FILE="workflow/Snakefile"
CONFIG_FILE="$DEFAULT_CONFIG_FILE"
RUN_DIR=""
SM_PARAMETERS=""
DRYRUN_FLAG=""
PROFILE=""


# ======================================================================================================================
usage="$(basename "$0") [-h] [-n] [-s SNAKEFILE] [-c CONFIG_FILE] [-w WORKDIR] [-P PROFILE] [-p \"snakemake parameters\"]
Launch the RIDE Snakemake workflow.

Options:
    -h                Show this help message and exit
    -n                Dry-run mode (equivalent to --dry-run)
    -s SNAKEFILE      Optional Snakefile (default: ${SNAKE_FILE})
    -c CONFIG_FILE    Config file to use (default: ${CONFIG_FILE})
    -w WORKDIR        Run directory:
                        • empty → ${RUNS_DIR}/<timestamp>/
                        • name  → ${RUNS_DIR}/<name>/
                        • path  → used as-is
    -P PROFILE        Snakemake profile (e.g. drmaa)
    -p PARAMETERS     Additional Snakemake parameters
"
# ======================================================================================================================

# ======================
# Parse CLI arguments
# ======================
while getopts ':hnc:p:w:s:P:' option; do
  case "$option" in
    h) echo "$usage"; exit 0 ;;
    n) DRYRUN_FLAG="--dry-run" ;;
    s) SNAKE_FILE=$OPTARG ;;
    c) CONFIG_FILE=$OPTARG ;;     # now strictly a file
    w) RUN_DIR=$OPTARG ;;
    P) PROFILE=$OPTARG ;;
    p) SM_PARAMETERS=$OPTARG ;;
    :) printf "Missing argument for -%s\n" "$OPTARG" >&2; echo "$usage" >&2; exit 1 ;;
    \?) printf "Illegal option: -%s\n" "$OPTARG" >&2; echo "$usage" >&2; exit 1 ;;
  esac
done
shift $((OPTIND - 1))

PROFILE_FLAG=""
if [[ -n "$PROFILE" ]]; then
    PROFILE_FLAG="--profile $PROFILE"
fi


# ===============================================
# Determine run directory logic (using RUNS_DIR)
# ===============================================
if [[ -z "$RUN_DIR" ]]; then
    TS=$(date +%F_%H%M%S)
    RUN_DIR="${RUNS_DIR}/${TS}"
else
    if [[ "$RUN_DIR" == */* ]]; then
        RUN_DIR="$RUN_DIR"           # full or relative path provided
    else
        RUN_DIR="${RUNS_DIR}/${RUN_DIR}"
    fi
fi

mkdir -p "$RUN_DIR"
RUN_DIR_PATH=$(readlink -f "$RUN_DIR")


# ====================
# Validate Snakefile
# ====================
if [[ ! -f "$SNAKE_FILE" ]]; then
    echo "ERROR: Snakefile not found: $SNAKE_FILE"
    exit 1
fi
SNAKEFILE_PATH=$(readlink -f "$SNAKE_FILE")


# ==============================================
# Validate config file (must be a single file)
# ==============================================
if [[ ! -f "$CONFIG_FILE" ]]; then
    echo "ERROR: Config file not found: $CONFIG_FILE"
    exit 1
fi
CONFIG_FILEPATH=$(readlink -f "$CONFIG_FILE")


# =============================================
# Copy full config folder for reproducibility
# =============================================
RUN_CONFIG_DIR="$RUN_DIR_PATH/config"
mkdir -p "$RUN_CONFIG_DIR"

if [[ -d "config" ]]; then
    cp -r config/* "$RUN_CONFIG_DIR/"
fi

# And ensure the used config file is copied explicitly over
cp "$CONFIG_FILEPATH" "$RUN_CONFIG_DIR/config.yaml"
CONFIG_USED="$RUN_CONFIG_DIR/config.yaml"


# ===========================================
# Activate environment (supports env paths)
# ===========================================
ENV_NAME="ride"
ENV_METADATA=".ride_envpath"

eval "$(conda shell.bash hook)"

if [[ -f "$ENV_METADATA" ]]; then
    ENVPATH=$(cat "$ENV_METADATA")

    if [[ -d "$ENVPATH" ]]; then
        echo "Activating environment via path: $ENVPATH"
        conda activate "$ENVPATH"
    else
        echo "WARNING: Env path not found, falling back to env name: $ENV_NAME"
        conda activate "$ENV_NAME"
    fi
else
    # auto-detect
    ENVPATH_DETECTED=$(conda env list | awk -v n="$ENV_NAME" '$1==n{print $NF}')
    if [[ -n "$ENVPATH_DETECTED" ]]; then
        echo "Activating environment detected at: $ENVPATH_DETECTED"
        conda activate "$ENVPATH_DETECTED"
    else
        echo "Activating environment by name: $ENV_NAME"
        conda activate "$ENV_NAME"
    fi
fi


# ------------------------------------------------------------
# Summary
# ------------------------------------------------------------
echo ""
echo "==============================================="
echo "            RIDE Snakemake launcher"
echo "==============================================="
echo "Environment     : $ENV_NAME"
echo "Run directory   : $RUN_DIR_PATH"
echo "Config file     : $CONFIG_USED"
echo "Snakefile       : $SNAKEFILE_PATH"

if [[ -n "$PROFILE" ]]; then
    echo "Profile         : $PROFILE"
else
    echo "Profile         : none"
fi

echo "Dry-run         : ${DRYRUN_FLAG:-no}"

if [[ -n "$SM_PARAMETERS" ]]; then
    echo "Extra params    : $SM_PARAMETERS"
else
    echo "Extra params    : none"
fi

echo "==============================================="
echo ""

# ------------------------------------------------------------
# Snakemake execution
# ------------------------------------------------------------
snakemake \
    --snakefile "$SNAKEFILE_PATH" \
    --directory "$RUN_DIR_PATH" \
    --use-conda \
    --configfile "$CONFIG_USED" \
    --summary  \
    --report "$RUN_DIR_PATH/report.html" \
    --printshellcmds \
    --restart-times 3 \
    --keep-going \
    --rerun-incomplete \
    $PROFILE_FLAG \
    --jobname "${ENV_NAME}.{rulename}.{jobid}.sh" \
    $DRYRUN_FLAG \
    $SM_PARAMETERS


echo ""
echo "RIDE pipeline completed."
echo ""
