#######################
import errno
import pandas as pd
import os
import multiprocessing
import psutil
from snakemake.utils import validate

report: "../report/workflow.rst"

validate(config, schema="../schemas/config.schema.yaml")

samples = pd.read_csv(config.get("samples"), sep='\t')
units = pd.read_csv(config.get("units"), sep='\t')
reheader = pd.read_csv(config.get("reheader"), sep='\t')

units_se=units[units["fq2"].isna()]
units_pe=units[units["fq2"].notna()]

def resolve_single_filepath(basepath, filename):
    return os.path.join(basepath, filename)

def resolve_results_filepath(dirname, filename):
    path = os.path.join(config.get('paths').get('results_dir'), dirname)
    return resolve_single_filepath(path, filename)

def resolve_logs_filepath(dirname, filename):
    path = os.path.join(config.get('paths').get('results_dir'), 'logs', dirname)
    return resolve_single_filepath(path, filename)

def resolve_benchmarks_filepath(dirname, filename):
    path = os.path.join(config.get('paths').get('results_dir'), 'benchmarks', dirname)
    return resolve_single_filepath(path, filename)

def resolve_envs_filepath(filename):
    path = os.path.join(config.get('paths').get('workdir'), 'workflow', 'envs')
    return resolve_single_filepath(path, filename)

def resolve_scripts_filepath(filename):
    path = os.path.join(config.get('paths').get('workdir'), 'workflow', 'scripts')
    return resolve_single_filepath(path, filename)

def ref_path(section, field):
    """
    Return the absolute path to a reference file defined in config.yaml.
    Path structure:
        <basepath>/<provider>/<release>/<file>
    """
    info = config["resources"][section]
    base = info["basepath"]
    provider = info.get("provider", "")
    release = info.get("release", "")
    filename = info[field]

    # Build structured path: base/provider/release/filename
    refdir = os.path.join(base, provider, release)

    return os.path.join(refdir, filename)

def temp_path(path=None):
    """
    Return a valid temporary directory path.
    """

    # Determine default temporary directory
    results_dir = config["paths"]["results_dir"]
    default_path = os.path.join(results_dir, "tmp")

    # If no explicit path is provided, return the default temp directory
    if not path:
        os.makedirs(default_path, exist_ok=True)
        return os.path.abspath(default_path)

    # Use custom path
    try:
        os.makedirs(path, exist_ok=True)
        return os.path.abspath(path)
    except Exception:
        # Fallback to default
        os.makedirs(default_path, exist_ok=True)
        return os.path.abspath(default_path)

def expand_filepath(filepath):
    filepath = os.path.expandvars(os.path.expanduser(filepath))
    if not os.path.isabs(filepath):
        raise FileNotFoundError(
            errno.ENOENT,
            os.strerror(errno.ENOENT) + " (path must be absolute)",
            filepath,
        )
    return filepath

def get_fastq(wildcards, units):
    if units.loc[wildcards.unit, ["fq2"]].isna().all():
        return expand_filepath(units.loc[wildcards.unit, ["fq1"]].dropna()[0])
    else:
        return expand_filepath(
            units.loc[wildcards.unit, ["fq1"]].dropna()[0]
        ), expand_filepath(units.loc[wildcards.unit, ["fq2"]].dropna()[0])


def get_a_fastq(wildcards, units, fq="fq1"):
    return expand_filepath(units.loc[wildcards.unit, [fq]].dropna()[0])

def get_sample_units(sample, layout="pe"):
    """
    Return the list of unit IDs (index values) for a given sample and layout.
    layout: "pe" or "se".
    """
    df = units_pe if layout == "pe" else units_se
    return df[df["sample"] == sample].index.tolist()

def get_unit_fastqs_pe(wildcards, read_pair="fq1"):
    """
    Return list of PE fastq files (fq1 or fq2) for a given sample.
    """
    unit_ids = get_sample_units(wildcards.sample, layout="pe")
    fastqs = []
    for u in unit_ids:
        val = units.loc[u, read_pair]
        if pd.notna(val):
            fastqs.append(expand_filepath(val))
    return fastqs


def get_unit_fastqs_se(wildcards, read_pair="fq1"):
    """
    Return list of SE fastq files (only fq1) for a given sample.
    """
    unit_ids = get_sample_units(wildcards.sample, layout="se")
    fastqs = []
    for u in unit_ids:
        val = units.loc[u, read_pair]
        if pd.notna(val):
            fastqs.append(expand_filepath(val))
    return fastqs


def get_unit_fastqs(wildcards, samples_df, label="units", read_pair="fq1"):
    """
    Return a list of FASTQ files for a given sample and a given read pair (fq1 or fq2).
    Units that do not contain the requested read (e.g., SE units when requesting fq2)
    are automatically skipped.
    """

    # Retrieve the comma-separated list of unit IDs for this sample
    units_str = samples_df.loc[wildcards.sample, label]
    unit_ids = units_str.split(",")

    fastqs = []
    for u in unit_ids:
        # Retrieve the FASTQ path for the requested read (fq1 or fq2)
        val = units.loc[u, read_pair]

        # Skip units that do not have this read (e.g., SE entries when read_pair == fq2)
        if pd.notna(val):
            fastqs.append(expand_filepath(val))

    return fastqs


def cpu_count():
    return multiprocessing.cpu_count()

def conservative_cpu_count(reserve_cores=1, max_cores=8):
    cores = max_cores if cpu_count() > max_cores else cpu_count()
    return max(cores - reserve_cores, 1)
