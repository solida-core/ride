import multiprocessing
import os.path
import os
import errno

def cpu_count():
    return multiprocessing.cpu_count()


def pipeline_cpu_count(reserve_cpu=2):
    cpu_nums = cpu_count() if reserve_cpu > cpu_count() \
        else cpu_count() - reserve_cpu
    return cpu_nums


def conservative_cpu_count(reserve_cores=1, max_cores=5):
    cores = max_cores if cpu_count() > max_cores else cpu_count()
    return max(cores - reserve_cores, 1)

def expand_filepath(filepath):
    filepath = os.path.expandvars(os.path.expanduser(filepath))
    if not os.path.isabs(filepath):
        raise FileNotFoundError(
            errno.ENOENT, os.strerror(errno.ENOENT)+" (path must be absolute)", filepath)
    return filepath

def references_abs_path(ref='references'):
    references = config.get(ref)
    basepath = expand_filepath(references['basepath'])
    provider = references['provider']
    release = references['release']

    return [os.path.join(basepath, provider, release)]



def resolve_single_filepath(basepath, filename):
    return [os.path.join(basepath, filename)]

def resolve_multi_filepath(basepath, dictionary):
    for k, v in dictionary.items():
        dictionary[k] = os.path.join(basepath, v)
    return dictionary


def get_references_label(ref='references'):
    references = config.get(ref)
    provider = references['provider']
    genome = references['release']

    return '_'.join([provider, genome])

def get_fastq(wildcards,samples,read_pair='fq'):
    return samples.loc[wildcards.sample,
                     [read_pair]].dropna()[0]

def get_cwd():
    return os.getcwd()