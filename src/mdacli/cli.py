#!/usr/bin/env python3
# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
#
# Copyright (c) 2021 Authors and contributors
#
# Released under the GNU Public Licence, v2 or any higher version
# SPDX-License-Identifier: GPL-2.0-or-later

"""
Main entry point for the MDAnalysis CLI interface.

This also demonstrates how other third party libraries could incorporate
this functionality.
"""
import argparse
import logging
import os
import sys
import traceback
import warnings

import MDAnalysis as mda
from MDAnalysis.analysis import __all__

from mdacli import __version__
from mdacli.colors import Emphasise
from mdacli.libcli import (
    KwargsDict,
    add_cli_universe,
    add_output_group,
    add_run_group,
    find_AnalysisBase_members_ignore_warnings,
    split_argparse_into_groups,
    )
from mdacli.logger import setup_logging
from mdacli.save import save_results
from mdacli.utils import convert_str_time, parse_callable_signature, parse_docs


logger = logging.getLogger(__name__)

# modules in MDAnalysis.analysis packages that are ignored by mdacli
# relevant modules used in this CLI factory
# hydro* are removed here because they have a different folder/file structure
# and need to be investigated separately
_skip_mods = ['base', 'hydrogenbonds', 'hbonds']

# skip modules that contain lists as parameters for AtomGroups since
# we can not parse them correctly so far.
_skip_mods += ['Contacts', 'Dihedral', 'PersistenceLength', 'InterRDF_s']

_relevant_modules = (mod for mod in __all__ if mod not in _skip_mods)

# serves CLI factory
STR_TYPE_DICT = {
    "bool": bool,
    "str": str,
    "list": list,
    "tuple": tuple,
    "dict": dict,
    "int": int,
    "float": float,
    "complex": complex,
    "NoneType": type(None),
    "AtomGroup": mda.AtomGroup,
    "Universe": mda.Universe,
    }


def _warning(message, *args, **kwargs):
    logger.warning(Emphasise.warning(message))


warnings.showwarning = _warning


def create_CLI(cli_parser, interface_name, parameters):
    """
    Add subparsers to `cli_parser`.

    Subparsers parameters are divided in the following categories:

    2. Analysis Run parameters
         time frame as begin, end, step and vebosity

    3. Saving Parameters
        output_prefix and output_directory

    4. Mandatory Parameters
        mandatory parameters are defined in the CLI as named parameters
        as per design

    5. Optional Parameters
        Named parameters in the Analysis class

    6. Reference Universe Parameters
        A reference Universe for selection commands. Only is created if
        AtomGroup arguments exist.

    All CLI's parameters are named parameters.

    Parameters
    ----------
    cli_parser : argparse.sub_parser
        The main parser where the new parser will be added.

    interface_name : str
        Name of the interface name.

    parameters : dict
        Parameters needed to fill the argparse requirements for the
        CLI interface.

    Returns
    -------
    None
    """
    analysis_class_parser = add_parser(cli_parser, interface_name, parameters)
    add_run_group(analysis_class_parser)

    # adds only if `save` method does not exist
    if not getattr(parameters['callable'], 'save', False):
        add_output_group(analysis_class_parser)

    # add positional and optional arguments
    pos_ = sorted(list(parameters["positional"].items()), key=lambda x: x[0])
    opt_ = sorted(list(parameters["optional"].items()), key=lambda x: x[0])

    parameters_to_parse = pos_ + opt_
    mandatory_parameters_group = analysis_class_parser.add_argument_group(
        title="Mandatory Parameters",
        description="Mandatory parameters of this Analysis",
        )
    groups = len(pos_) * [mandatory_parameters_group]

    # Only create parser if optional arguments exist
    if len(opt_) > 0:
        optional_parameters_group = analysis_class_parser.add_argument_group(
            title="Optional Parameters",
            description="Optional parameters specific of this Analysis",
            )
        groups += len(opt_) * [optional_parameters_group]

    for group, (name, args_dict) in zip(groups, parameters_to_parse):

        # prepares parameters before add_argument
        try:
            type_ = STR_TYPE_DICT[args_dict["type"]]
        except KeyError:
            type_ = str

        try:
            default = args_dict["default"]
        except KeyError:
            default = None

        name_par = "-" + name
        description = args_dict["desc"]
        if issubclass(type_, (list, tuple)):
            group.add_argument(
                name_par,
                dest=name,
                default=default,
                nargs="+",
                help="{} (default: %(default)s)".format(description)
                )
        elif issubclass(type_, dict):
            group.add_argument(
                name_par,
                dest=name,
                default=None,
                action=KwargsDict,
                help=description,
                )
        elif type_ is bool:
            group.add_argument(
                name_par,
                dest=name,
                action="store_false" if default else "store_true",
                default=default,
                help=description,
                )
        elif type_ is mda.AtomGroup:
            group.add_argument(
                name_par,
                dest=name,
                type=str,
                default=default,
                help=description + " Use a MDAnalysis selection string."
                )

            # Create one reference Universe argument for atom selection
            try:
                reference_universe_group
            except NameError:
                reference_universe_group = \
                    analysis_class_parser.add_argument_group(
                        title="Reference Universe Parameters",
                        description="Parameters specific for loading "
                                    "the reference topology and trajectory"
                                    " used for atom selection.")
                add_cli_universe(reference_universe_group)

        elif issubclass(type_, mda.Universe):
            add_cli_universe(group, name)
        else:
            group.add_argument(
                name_par,
                dest=name,
                type=type_,
                default=default,
                help=f"{description} (default: %(default)s)"
                )
    return


def add_parser(cli_parser, interface_name, parameters):
    """Add the parser.

    str, dict -> None
    ...
    """
    # creates the subparser
    analysis_class_parser = cli_parser.add_parser(
        interface_name,
        help=parameters["desc"],
        description=f"{parameters['desc']}\n\n{parameters['desc_long']}",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        )

    # Add run_analsis function as the default func parameter.
    # this is possible because the run_analsis function is equal to all
    # Analysis Classes
    analysis_class_parser.set_defaults(
        analysis_callable=parameters["callable"])

    return analysis_class_parser


def create_universe(topology,
                    coordinates=None,
                    topology_format=None,
                    atom_style=None,
                    trajectory_format=None):
    """
    Initilize a MDAnalysis universe instance.

    Parameters
    ----------
    topology : str, stream, `~MDAnalysis.core.topology.Topology`, `np.ndarray`
        A CHARMM/XPLOR PSF topology file, PDB file or Gromacs GRO file; used to
        define the list of atoms. If the file includes bond information,
        partial charges, atom masses, ... then these data will be available to
        MDAnalysis. Alternatively, an existing
        :class:`MDAnalysis.core.topology.Topology` instance may be given,
        numpy coordinates, or None for an empty universe.
    coordinates : str, stream, list of str, list of stream
        Coordinates can be provided as files of
        a single frame (eg a PDB, CRD, or GRO file); a list of single
        frames; or a trajectory file (in CHARMM/NAMD/LAMMPS DCD, Gromacs
        XTC/TRR, or generic XYZ format). The coordinates must be
        ordered in the same way as the list of atoms in the topology.
        See :ref:`Supported coordinate formats` for what can be read
        as coordinates. Alternatively, streams can be given.
    topology_format : str, None
        Provide the file format of the topology file; ``None`` guesses it from
        the file extension. Can also pass a subclass of
        :class:`MDAnalysis.topology.base.TopologyReaderBase` to define a custom
        reader to be used on the topology file.
    trajectory_format : str or list or object
        provide the file format of the coordinate or trajectory file;
        ``None`` guesses it from the file extension. Note that this
        keyword has no effect if a list of file names is supplied because
        the "chained" reader has to guess the file format for each
        individual list member [``None``]. Can also pass a subclass of
        :class:`MDAnalysis.coordinates.base.ProtoReader` to define a custom
        reader to be used on the trajectory file.
    atom_style : str
        Customised LAMMPS `atom_style` information. Only works with
        `topology_format = data`

    Returns
    -------
    `MDAnalysis.Universe`
    """
    universe = mda.Universe(topology,
                            topology_format=topology_format,
                            atom_style=atom_style)

    if coordinates is not None:
        universe.load_new(coordinates, format=trajectory_format)

    return universe


def run_analsis(analysis_callable,
                mandatory_analysis_parameters,
                optional_analysis_parameters=None,
                reference_universe_parameters=None,
                run_parameters=None,
                output_parameters=None):
    """Perform main client logic.

    Parameters
    ----------
    analysis_callable : function
        Analysis class for which the analysis is performed.
    mandatory_analysis_parameters : dict
        Mandatory parameters for executing the analysis
    optional_analysis_parameters : dict
        Optional parameters for executing the analysis
    run_parameters : dict
        time frame parameters: start, stop, step, verbose
    output_parameters : dict
        output_prefix and output_directory

    Returns
    -------
    `MDAnalysis.analysis.base.AnalysisBase`
        AnalysisBase instance of the given ``analysis_callable`` after run.
    """
    if optional_analysis_parameters is None:
        optional_analysis_parameters = {}

    if run_parameters is None:
        run_parameters = {}

    if output_parameters is None:
        output_parameters = {}

    verbose = run_parameters.pop("verbose", False)

    if reference_universe_parameters is not None:
        reference_universe = create_universe(**reference_universe_parameters)
    else:
        reference_universe = None

    # Initilize analysis callable
    universe = convert_analysis_parameters(analysis_callable,
                                           mandatory_analysis_parameters,
                                           reference_universe)

    convert_analysis_parameters(analysis_callable,
                                optional_analysis_parameters,
                                reference_universe)

    if universe is None:
        universe = reference_universe

    ac = analysis_callable(**mandatory_analysis_parameters,
                           **optional_analysis_parameters)

    # Run the analysis
    for key, value in run_parameters.items():
        run_parameters[key] = convert_str_time(value, universe.trajectory.dt)

    ac.run(verbose=verbose, **run_parameters)

    # Save results
    try:
        ac.save_results()

    except AttributeError:
        directory = output_parameters.get("output_directory", "")
        fname = output_parameters.get("output_prefix", "")
        fname = f"{fname}_{analysis_callable.__name__}" if fname \
            else analysis_callable.__name__
        save_results(ac.results, os.path.join(directory, fname))

    return ac


def convert_analysis_parameters(analysis_callable,
                                analysis_parameters,
                                reference_universe=None):
    """
    Convert parameters from the command line suitable for anlysis.

    Special types (i.e AtomGroups, Universes) are converted from the command
    line strings into the correct format. Parameters are changed inplace.
    Note that only keys are converted and no new key are added if
    present in the doc of the `analysis_callable` but not
    in the `analysis_parameters` dict.

    The following types are converted:

    * AtomGroup: Select atoms based on ``universe.select_atoms``
    * Universe: Created from parameters.

    Parameters
    ----------
    analysis_callable : function
        Analysis class for which the analysis should be performed.
    analysis_parameters : dict
        parameters to be processed
    reference_universe : `MDAnalysis.Universe`
        Universe from which the AtomGroup selection are done.

    Returns
    -------
    universe : Universe
        The universe created from the anaylysis parameters or None
        of no ine is created

    Raises
    ------
    ValueError
        If an Atomgroup does not contain any atoms
    """
    params = parse_docs(analysis_callable)[2]
    universe = None

    # If a Universe is part of the parameters several extra arguments with
    # non matching names were created. We seperate them by their connecting
    # character.
    analysis_parameters_keys = [p.split("-")[-1] for p
                                in analysis_parameters.keys()]

    for param_name, dictionary in params.items():
        if param_name in analysis_parameters_keys:
            if "AtomGroup" in dictionary['type']:
                sel = reference_universe.select_atoms(
                    analysis_parameters[param_name])
                if sel:
                    analysis_parameters[param_name] = sel
                else:
                    raise ValueError(f"AtomGroup `-{param_name}`"
                                     f" with selection"
                                     f" `{analysis_parameters[param_name]}`"
                                     f" does not contain any atoms")
            elif "Universe" in dictionary['type']:

                # All Universe parameters
                universe_parameters = {"topology": None,
                                       "coordinates": None,
                                       "topology_format": None,
                                       "atom_style": None,
                                       "trajectory_format": None}

                for k in universe_parameters.keys():
                    universe_parameters[k] = analysis_parameters.pop(
                        f"{k}-{param_name}")

                universe = create_universe(**universe_parameters)
                analysis_parameters[param_name] = universe

    return universe


def setup_clients(ap, title, members):
    """
    Set up analysis clients for an ArgumentParser instance.

    Parameters
    ----------
    ap : `argparse.ArgumentParser`
        Argument parser instance
    title : str
        title of the parser
    members : list
        list containing Analysis classes for setting up the parser
    """
    cli_subparser = ap.add_subparsers(title=title)

    analysis_interfaces = {
        member.__name__: parse_callable_signature(member)
        for member in members
        }

    # adds each Analysis class/function as a CLI under 'cli_subparser'
    # to be writen
    for interface_name, parameters in analysis_interfaces.items():
        create_CLI(cli_subparser, interface_name, parameters)


def init_base_argparse():
    """Create basic argument parser."""
    ap = argparse.ArgumentParser()

    ap.add_argument(
        '--version',
        action='version',
        version="mdacli {}".format(__version__),
        )

    ap.add_argument(
        '--debug',
        action='store_true',
        help="Run with debug options.",
        )

    ap.add_argument('--logfile',
                    dest='logfile',
                    action='store',
                    help='Logfile (optional)')
    return ap


def _exit_if_a_is_b(obj1, obj2, msg):
    """Exit if `obj1` and `obj2` are the same."""
    if obj1 is obj2:
        sys.exit(msg)


def main():
    """Execute main CLI entry point."""
    modules = find_AnalysisBase_members_ignore_warnings(_relevant_modules)
    _exit_if_a_is_b(modules, None, "No analysis modules founds.")

    ap = init_base_argparse()

    if len(sys.argv) < 2:
        ap.error("A subcommand is required.")

    # There is to much useless code execution done here:
    # 1. We do not have to setup all possible clients all the time.
    #    i.e. for `mdacli RMSD` only the RMSD client should be build.
    # 2. for something like `mdacli -h` We do not have to build every
    #   sub parser in complete detail.
    setup_clients(ap, title="MDAnalysis Analysis Clients", members=modules)

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
            _key = analysis_callable.__name__
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
