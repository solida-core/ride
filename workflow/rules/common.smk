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

def temp_path(path=""):
    default_path = os.path.join(config.get('paths').get('results_dir'), 'tmp')
    if path:
        try:
            os.makedirs(path)
        except OSError as e:
            if e.errno != errno.EEXIST:
                return default_path
        return path
    return default_path

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

def get_unit_fastqs(wildcards, samples, label='units',read_pair='fq'):
    for unit_set in samples.loc[wildcards.sample,[label]]:
        print(wildcards.sample)
    return [expand_filepath(units.loc[x,[read_pair]].dropna()[0]) for x in unit_set.split(',')]


def cpu_count():
    return multiprocessing.cpu_count()

def conservative_cpu_count(reserve_cores=1, max_cores=8):
    cores = max_cores if cpu_count() > max_cores else cpu_count()
    return max(cores - reserve_cores, 1)
