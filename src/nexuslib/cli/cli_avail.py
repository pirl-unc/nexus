# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.


"""
The purpose of this python3 script is to create parser
and run nexus 'avail' command.
"""


import argparse
import glob
import math
import os
import random
from importlib import resources
from pathlib import Path
from collections import OrderedDict


def add_nexus_avail_arg_parser(sub_parsers):
    """
    Adds 'avail' parser.

    Parameters
    ----------
    sub_parsers     :   argparse.ArgumentParser subparsers object.

    Returns
    -------
    sub_parsers     :   argparse.ArgumentParser subparsers object.
    """
    parser = sub_parsers.add_parser(
        'avail',
        help='Lists all available workflows.'
    )
    parser._action_groups.pop()
    parser.set_defaults(which='avail')
    return sub_parsers


def run_nexus_avail_from_parsed_args(args):
    """
    Runs Nexus 'avail' command using parameters from parsed arguments.
    """
    print("Available workflows:")
    resources_path = resources.files('nexuslib').joinpath("pipelines")
    workflows_dict = OrderedDict()
    nf_scripts = glob.glob("%s/**/*.nf" % resources_path, recursive=True)
    for nf_script in nf_scripts:
        nf_script_dir = nf_script.replace(str(resources_path), '')
        nf_script_dirs = nf_script_dir.split('/')
        if nf_script_dirs[1] == 'modules':
            continue
        if nf_script_dirs[1] not in workflows_dict.keys():
            workflows_dict[nf_script_dirs[1]] = []
        workflows_dict[nf_script_dirs[1]].append(os.path.basename(nf_script))
    last_key = ''
    for key, val in workflows_dict.items():
        if last_key != key:
            last_key = key
            print('\t%s:' % (key.replace('_', ' ')[0].capitalize() + key.replace('_', ' ')[1:]))
        for v in val:
            print('\t\t%s' % v)
