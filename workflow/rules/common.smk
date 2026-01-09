import os
import sys
import errno
import pandas as pd
import multiprocessing
import psutil
from snakemake.utils import validate


###############################################
# Validate config schema
###############################################
report: "../report/workflow.rst"
validate(config, schema="../schemas/config.schema.yaml")

###############################################
# Load input tables
###############################################

samples = pd.read_table(
    config["samples"],
    sep=r"\s+",
    dtype=str,
    keep_default_na=True,
    na_values=["", " ", "NA", "NaN", "nan", "NONE", "None"]
).set_index("sample", drop=False)

validate(samples.to_dict(orient="list"), schema="../schemas/samples.schema.yaml")


units = pd.read_table(
    config["units"],
    sep=r"\s+",
    dtype=str,
    keep_default_na=True,
    na_values=["", " ", "NA", "NaN", "nan", "NONE", "None"]
).set_index("unit", drop=False)

units = units.where(pd.notna(units), None)

validate(units.to_dict(orient="list"), schema="../schemas/units.schema.yaml")


try:
    reheader = pd.read_table(
        config["reheader"],
        sep=r"\s+",
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
# Path resolvers
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
# Reference resolver
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
# Temporary path handling
###############################################
def temp_path(path=None):
    default_path = os.path.join("tmp")

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
# FASTQ helpers
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

