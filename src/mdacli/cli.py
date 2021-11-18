#!/usr/bin/env python3
# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
#
# Copyright (c) 2021 Authors and contributors
#
# Released under the GNU Public Licence, v2 or any higher version
# SPDX-License-Identifier: GPL-2.0-or-later

"""The toplevel command line interface."""
import logging
import sys
import traceback
import warnings

from .colors import Emphasise
from .libcli import (
    find_AnalysisBase_members,
    init_base_argparse,
    run_analsis,
    setup_clients,
    split_argparse_into_groups,
    )
from .logger import setup_logging
from .utils import _exit_if_a_is_b


logger = logging.getLogger(__name__)


def cli(name,
        module_list,
        version="",
        description="",
        skip_modules=None,
        ignore_warnings=False):
    """Create the command-line interface.

    This function creates a command line interface with a given `name` based
    on a module list.

    Parameters
    ----------
    name : str
        name of the interface
    module_list : list
        list of module from which the cli is build up.
    description : str
        description of the cli
    skip_modules : list
        list of strings containing modules that should be ommited
    ignore_warnings : bool
        ignore warnings when importing modules

    Examples
    --------
    This example code creates the command line interface of MDAnalysis::

        from MDAnalysis.analysis import __all__

        import mdacli

        skip_mods = ['AnalysisFromFunction',
                    'HydrogenBondAnalysis',
                    'WaterBridgeAnalysis',
                    'Contacts',
                    'Dihedral',
                    'PersistenceLength',
                    'InterRDF_s']

        mdacli.cli(name="MDAnalysis",
                   module_list=__all__,
                   version=mdacli.__version__,
                   description=__doc__,
                   skip_modules=skip_mods,
                   ignore_warnings=True)
    """
    modules = find_AnalysisBase_members(module_list,
                                        ignore_warnings=ignore_warnings)

    modules = [mod for mod in modules if mod.__name__ not in skip_modules]
    _exit_if_a_is_b(modules, None, "No analysis modules founds.")

    ap = init_base_argparse(name=name,
                            version=version,
                            description=description)

    if len(sys.argv) < 2:
        ap.error("A subcommand is required.")

    # There is to much useless code execution done here:
    # 1. We do not have to setup all possible clients all the time.
    #    i.e. for `mda RMSD` only the RMSD client should be build.
    # 2. for something like `mdacli -h` We do not have to build every
    #   sub parser in complete detail.
    setup_clients(ap, title=f"{name} Analysis Modules", members=modules)

    # Be case insensitive for the subcommand
    sys.argv[1] = sys.argv[1].lower()

    args = ap.parse_args()

    if args.debug:
        args.verbose = True
    else:
        # Ignore all warnings if not in debug mode
        warnings.filterwarnings("ignore")

    with setup_logging(logger, logfile=args.logfile, debug=args.debug):
        # Execute the main client interface.
        try:
            analysis_callable = args.analysis_callable

            # Get the correct ArgumentParser instance from all subparsers
            # `[0]` selects the first subparser where our analysises live in.
            _key = analysis_callable.__name__.lower()
            ap_sup = ap._subparsers._group_actions[0].choices[_key]
            arg_grouped_dict = split_argparse_into_groups(ap_sup, args)

            # Some parameters may not exist
            arg_grouped_dict.setdefault("Optional Parameters", {})
            arg_grouped_dict.setdefault("Reference Universe Parameters", None)

            run_analsis(analysis_callable,
                        arg_grouped_dict["Mandatory Parameters"],
                        arg_grouped_dict["Optional Parameters"],
                        arg_grouped_dict["Reference Universe Parameters"],
                        arg_grouped_dict["Analysis Run Parameters"],
                        arg_grouped_dict["Output Parameters"])
        except Exception as e:
            if args.debug:
                traceback.print_exc()
            else:
                sys.exit(Emphasise.error(f"Error: {e}"))
