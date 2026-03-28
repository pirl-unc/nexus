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
from typing import List


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


def get_available_workflows():
    """
    Fetches all available workflows. Scans both 'subworkflows' (single-tool workflows) and 'workflows'.
    Process definitions in 'tools' are excluded.

    Returns
    -------
    nf_scripts_paths    :   Dictionary (key = workflow name, value = path to nextflow script).
    """
    nf_scripts_paths = {}
    base_path = resources.files('nexuslib')
    for directory in ('subworkflows', 'workflows'):
        search_path = base_path.joinpath(directory)
        nf_scripts = glob.glob("%s/**/*.nf" % search_path, recursive=True)
        for nf_script in nf_scripts:
            workflow_name = os.path.basename(nf_script)
            nf_scripts_paths[workflow_name] = nf_script
    return nf_scripts_paths


def run_workflow(
        workflow: str,
        nextflow: str,
        workflow_args: List[str]
):
    """
    Runs a workflow.

    Args
    ----
    workflow        :   Workflow (e.g. 'long_read_alignment_minimap2').
    nextflow        :   Nextflow path.
    workflow_args   :   List of workflow specific arguments.
    """
    nf_scripts_paths = get_available_workflows()
    if workflow not in nf_scripts_paths:
        raise Exception(f"{workflow} is not available.")

    if nextflow == "nextflow":
        nextflow = get_alias_path(executable_name="nextflow")

    command = [nextflow, "run", nf_scripts_paths[workflow], "-resume"] + workflow_args

    print(' '.join(command), flush=True)

    process = subprocess.Popen(
        ' '.join(command),
        shell=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        bufsize=1,
        executable='/bin/bash'
    )

    for line in process.stdout:
        print(line, end="", flush=True)

    return_code = process.wait()

    if return_code == 0:
        print("Nextflow returned zero (ran successfully).", flush=True)
    else:
        raise Exception('Nextflow returned non-zero exit code: %i' % return_code)
