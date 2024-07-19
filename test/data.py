import os
import subprocess


def get_data_path(name: str) -> str:
    """
    Return the absolute path to a file in the test/data directory.

    Parameters
    ----------
    name    :   Name of file.

    Returns
    -------
    Absolute path to a file in the test/data directory.
    """
    return os.path.join(os.path.dirname(__file__), "data", name)


def get_alias_path(executable_name: str):
    try:
        result = subprocess.run('which %s' % executable_name,
                                stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE,
                                shell=True,
                                text=True,
                                executable='/bin/bash')
        return result.stdout.strip()
    except:
        raise Exception('%s not found.' % executable_name)


def get_path(path: str):
    try:
        result = subprocess.run('echo $%s' % path,
                                stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE,
                                shell=True,
                                text=True,
                                executable='/bin/bash')
        return result.stdout.strip()
    except:
        raise Exception('%s not found.' % path)
