###############################################
# 1) Imports
###############################################
import os
import sys
import errno
import pandas as pd
import multiprocessing
import psutil
from snakemake.utils import validate


###############################################
# 2) Validate config schema
###############################################
report: "../report/workflow.rst"
validate(config, schema="../schemas/config.schema.yaml")


###############################################
# 3) Safe TSV loading
###############################################
def safe_read_tsv(path, required_columns=None, allow_empty=False):
    """
    Safely read a TSV file with:
      - existence check
      - parsing check
      - optional empty handling
      - column validation

        If allow_empty == True, empty file returns None.
    """

    # Expand env variables and ~
    path = os.path.expandvars(os.path.expanduser(path))
    filename = os.path.basename(path)

    # Check file existence
    if not os.path.isfile(path):
        print(f"\n[ERROR] File not found: {path}\n", file=sys.stderr)
        sys.exit(1)

    # Try reading TSV
    try:
        df = pd.read_csv(
            path,
            sep="\t",
            dtype=str,
            keep_default_na=True,
            na_values=["", " ", "NA", "NaN", "nan", "NONE", "None"]
        )
    except EmptyDataError:
        # File exists but has no content or header
        if allow_empty:
            print(f"[WARNING] {filename} is empty or has no parsable content → returning None")
            return None
        else:
            print(f"\n[ERROR] {filename} is empty or has no columns: {path}\n", file=sys.stderr)
            sys.exit(1)
    except Exception as e:
        print(f"\n[ERROR] Could not read {filename}: {path}", file=sys.stderr)
        print(f"Details: {e}\n", file=sys.stderr)
        sys.exit(1)

    # Empty file handling
    if df.empty:
        if allow_empty:
            print(f"[WARNING] {filename} is empty → returning None")
            return None
        else:
            print(f"\n[ERROR] {filename} is empty: {path}\n", file=sys.stderr)
            sys.exit(1)

    # Require at least 2 columns
    if df.shape[1] < 2:
        print(f"\n[ERROR] {filename} contains fewer than 2 columns.\n", file=sys.stderr)
        sys.exit(1)

    # Validate required columns
    if required_columns:
        required_columns = set(required_columns)
        missing = required_columns - set(df.columns)
        if missing:
            print(f"\n[ERROR] {filename} is missing required columns:", file=sys.stderr)
            for col in missing:
                print(f"  - {col}", file=sys.stderr)
            sys.exit(1)

    print(f"[OK] Loaded {filename}: {df.shape[0]} rows, {df.shape[1]} columns")
    return df



###############################################
# Load input tables
###############################################

samples = pd.read_table(
    config["samples"],
    dtype=str,
    keep_default_na=True,
    na_values=["", " ", "NA", "NaN", "nan", "NONE", "None"]
).set_index("sample", drop=False)

validate(samples.to_dict(orient="list"), schema="../schemas/samples.schema.yaml")


units = pd.read_table(
    config["units"],
    dtype=str,
    keep_default_na=True,
    na_values=["", " ", "NA", "NaN", "nan", "NONE", "None"]
).set_index("unit", drop=False)

validate(units.to_dict(orient="list"), schema="../schemas/units.schema.yaml")


try:
    reheader = pd.read_table(
        config["reheader"],
        dtype=str
    )
    validate(reheader.to_dict(orient="list"), schema="../schemas/reheader.schema.yaml")
except pd.errors.EmptyDataError:
    reheader = None



###############################################
# 5) Split SE / PE and sample lists
###############################################
units["fq2"] = units["fq2"].replace({"": pd.NA, " ": pd.NA})

units_se = units[units["fq2"].isna()].copy()
units_pe = units[units["fq2"].notna()].copy()

SAMPLES     = samples["sample"].tolist()
SAMPLES_PE  = units_pe["sample"].unique().tolist()
SAMPLES_SE  = units_se["sample"].unique().tolist()


###############################################
# 6) Path resolvers
###############################################
def resolve_single_filepath(basepath, filename):
    return os.path.join(basepath, filename)

def resolve_results_filepath(dirname, filename):
    return os.path.join("results", dirname, filename)

def resolve_logs_filepath(dirname, filename):
    return os.path.join("logs", dirname, filename)

def resolve_benchmarks_filepath(dirname, filename):
    return os.path.join("benchmarks", dirname, filename)

def resolve_envs_filepath(filename):
    return os.path.join(workflow.basedir, "envs", filename)

def resolve_scripts_filepath(filename):
    return os.path.join(workflow.basedir, "scripts", filename)



###############################################
# 7) Reference resolver
###############################################
def ref_path(section, field):
    """
    Construct a structured path:
    <basepath>/<provider>/<release>/<filename>
    """
    info = config["resources"][section]
    base = info["basepath"]
    provider = info.get("provider", "")
    release  = info.get("release", "")
    filename = info[field]

    return os.path.join(base, provider, release, filename)


###############################################
# 8) Temporary path handling
###############################################
def temp_path(path=None):
    results_dir = config["paths"]["results_dir"]
    default_path = os.path.join(results_dir, "tmp")

    if path is None:
        os.makedirs(default_path, exist_ok=True)
        return os.path.abspath(default_path)

    try:
        os.makedirs(path, exist_ok=True)
        return os.path.abspath(path)
    except Exception:
        os.makedirs(default_path, exist_ok=True)
        return os.path.abspath(default_path)


###############################################
# 9) FASTQ helpers
###############################################
def expand_filepath(filepath):
    filepath = os.path.expandvars(os.path.expanduser(filepath))
    if not os.path.isabs(filepath):
        raise FileNotFoundError(
            errno.ENOENT,
            os.strerror(errno.ENOENT) + " (path must be absolute)",
            filepath
        )
    return filepath


def get_sample_units(sample, layout="pe"):
    df = units_pe if layout == "pe" else units_se
    return df[df["sample"] == sample].index.tolist()


def get_unit_fastqs_pe(wildcards, read_pair="fq1"):
    unit_ids = get_sample_units(wildcards.sample, "pe")
    return [
        expand_filepath(units.loc[u, read_pair])
        for u in unit_ids if pd.notna(units.loc[u, read_pair])
    ]


def get_unit_fastqs_se(wildcards, read_pair="fq1"):
    unit_ids = get_sample_units(wildcards.sample, "se")
    return [
        expand_filepath(units.loc[u, read_pair])
        for u in unit_ids if pd.notna(units.loc[u, read_pair])
    ]


###############################################
# 10) CPU helpers
###############################################
def cpu_count():
    return multiprocessing.cpu_count()

def conservative_cpu_count(reserve_cores=1, max_cores=8):
    cores = min(cpu_count(), max_cores)
    return max(cores - reserve_cores, 1)

