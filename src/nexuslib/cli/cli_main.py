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
The purpose of this python3 script is to implement the primary nexus command.
"""


import argparse
import nexuslib
from .cli_run import *
from .cli_avail import *


def init_arg_parser():
    """
    Initializes the input argument parser.

    Returns
    -------
    An instance of argparse.ArgumentParser
    An instance of argparse.ArgumentParser subparsers
    """
    arg_parser = argparse.ArgumentParser(
        description="Nexus: NEXtflow's Ultimate Streamliner."
    )
    arg_parser.add_argument(
        '--version', '-v',
        action='version',
        version='%(prog)s version ' + str(nexuslib.__version__)
    )
    sub_parsers = arg_parser.add_subparsers(help='nexus sub-commands.')
    return arg_parser, sub_parsers


def run():
    # Step 1. Initialize argument parser
    arg_parser, sub_parsers = init_arg_parser()
    sub_parsers = add_nexus_avail_arg_parser(sub_parsers=sub_parsers)       # avail
    sub_parsers = add_nexus_run_arg_parser(sub_parsers=sub_parsers)         # run
    args, unknown = arg_parser.parse_known_args()

    # Step 2. Execute function based on CLI arguments
    if args.which == 'avail':
        run_nexus_avail_from_parsed_args(args=args)
    elif args.which == 'run':
        run_nexus_run_from_parsed_args(args=args, workflow_args=unknown)
    else:
        raise Exception("Invalid command: %s" % args.which)