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
and run nexus 'run' command.
"""


import argparse
from ..main import run_workflow


def add_nexus_run_arg_parser(sub_parsers):
    """
    Adds 'generate' parser.

    Parameters
    ----------
    sub_parsers     :   argparse.ArgumentParser subparsers object.

    Returns
    -------
    sub_parsers     :   argparse.ArgumentParser subparsers object.
    """
    parser = sub_parsers.add_parser(
        'run',
        add_help=False
    )
    parser._action_groups.pop()

    # Required arguments
    parser_required = parser.add_argument_group('required arguments')
    parser_required.add_argument(
        "--nf-workflow",
        dest="nf_workflow",
        type=str,
        help="Nextflow workflow. To view the list of all available workflows, run 'nexus avail'."
    )
    # Optional arguments
    parser_optional = parser.add_argument_group('optional arguments')
    parser_optional.add_argument(
        "--nextflow",
        dest="nextflow",
        type=str,
        default='nextflow',
        required=False,
        help="nextflow path (default: nextflow)."
    )
    parser.set_defaults(which='run')
    return sub_parsers


def run_nexus_run_from_parsed_args(args, workflow_args):
    """
    Runs nexus 'run' command using parameters from parsed arguments.

    Parameters
    ----------
    args    :   argparse.ArgumentParser object
                with the following variables:
                workflow
    """
    # Print --help if requested
    if args.nf_workflow is None and ('--help' in workflow_args or '-h' in workflow_args):
        print("usage: nexus run [-h] [--nf-workflow NF_WORKFLOW] [--nextflow NEXTFLOW] [workflow-specific parameters]")
        print("")
        print("required arguments:")
        print("\t--nf-workflow NF_WORKFLOW. To view the list of all available workflows, run 'nexus avail'")
        print("")
        print("optional arguments:")
        print("\t--nextflow NEXTFLOW (default: nextflow).'")
        print("")
        print("view workflow-specific parameters:")
        print("\tnexus run --nf-workflow NF_WORKFLOW --help")
    else:
        run_workflow(workflow=args.nf_workflow,
                     nextflow=args.nextflow,
                     workflow_args=workflow_args)
