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


import glob
import os
import subprocess
from importlib import resources
from pathlib import Path
from typing import List


def get_available_workflows():
    """
    Fetches all available workflows.

    Returns
    -------
    nf_scripts_paths    :   Dictionary (key = workflow name, value = path to nextflow script).
    """
    resources_path = resources.path('nexuslib', "pipelines")
    nf_scripts_paths = {}
    nf_scripts = glob.glob("%s/**/*.nf" % resources_path, recursive=True)
    for nf_script in nf_scripts:
        nf_script_dir = nf_script.replace(str(resources_path), '')
        nf_script_dirs = nf_script_dir.split('/')
        if nf_script_dirs[1] == 'modules':
            continue
        workflow_name = os.path.basename(nf_script)
        nf_scripts_paths[workflow_name] = nf_script
    return nf_scripts_paths


def run_workflow(workflow: str, workflow_args: List[str]):
    """
    Runs a workflow.

    Args
    ----
    workflow        :   Workflow (e.g. 'long_read_alignment_minimap2').
    workflow_args   :   List of workflow specific arguments.
    """
    nf_scripts_paths = get_available_workflows()
    if workflow in nf_scripts_paths:
        command = ['nextflow', 'run', nf_scripts_paths[workflow]] + workflow_args
        print(' '.join(command))
        process = subprocess.Popen(
            ' '.join(command),
            shell=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True
        )
        while True:
            output = process.stdout.readline()
            if output == '' and process.poll() is not None:
                break
            if output:
                print(output.strip())
        return_code = process.wait()
        print("Return Code:", return_code)
    else:
        raise Exception('%s is not available.' % workflow)
