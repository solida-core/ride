import multiprocessing
import os.path


def cpu_count():
    return multiprocessing.cpu_count()


def pipeline_cpu_count(reserve_cpu=2):
    cpu_nums = cpu_count() if reserve_cpu > cpu_count() \
        else cpu_count() - reserve_cpu
    return cpu_nums


def references_abs_path(ref='references'):
    references = config.get(ref)
    basepath = references['basepath']
    provider = references['provider']
    release = references['release']

    return [os.path.join(basepath, provider, release)]


def resolve_single_filepath(basepath, filename):
    return [os.path.join(basepath, filename)]


def get_references_label(ref='references'):
    references = config.get(ref)
    provider = references['provider']
    genome = references['release']

    return '_'.join([provider, genome])
