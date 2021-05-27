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
import os
import sys
import warnings

import MDAnalysis as mda
from MDAnalysis.analysis import __all__

from mdacli import __version__
from mdacli.colors import Emphasise
from mdacli.libcli import KwargsDict, find_AnalysisBase_members_ignore_warnings
from mdacli.save import save_results
from mdacli.utils import convert_str_time, parse_callable_signature, parse_docs


# modules in MDAnalysis.analysis packages that are ignored by mdacli
# relevant modules used in this CLI factory
# hydro* are removed here because they have a different folder/file structure
# and need to be investigated separately
_skip_mods = ('base', 'hydrogenbonds', 'hbonds')
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
    }


def _warning(message, *args, **kwargs):
    print(Emphasise.warning(f"Warning: {message}"))


warnings.showwarning = _warning


def create_CLI(cli_parser, interface_name, parameters):
    """
    Add subparsers to `cli_parser`.

    Subparsers parameters are divided in three categories:

    1) Common Analysis classes Parameters
        Common to all generated CLIs, for example:
            * topology, trajectory, time frame

    2) Mandatory Parameters
        mandatory parameters are defined in the CLI as named parameters
        as per design

    3) Optional Parameters
        Named parameters in the Analysis class

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

    universe_group = analysis_class_parser.add_argument_group(
        title="Universe Parameters",
        description="Parameters specific for loading the topology and"
                    " trajectory"
        )

    universe_group.add_argument(
        "-s",
        dest="topology",
        type=str,
        default="topol.tpr",
        help="The topolgy file. "
        "The FORMATs {} are implemented in MDAnalysis."
        "".format(", ".join(mda._PARSERS.keys())),
        )

    universe_group.add_argument(
        "-top",
        dest="topology_format",
        type=str,
        default=None,
        help="Override automatic topology type detection. "
        "See topology for implemented formats.")

    universe_group.add_argument(
        "-atom_style",
        dest="atom_style",
        type=str,
        default=None,
        help="Manually set the atom_style information"
        "(currently only LAMMPS parser). E.g. atom_style='id type x y z'.")

    universe_group.add_argument(
        "-f",
        dest="coordinates",
        type=str,
        default=None,
        nargs="+",
        help="A single or multiple coordinate files. "
        "The FORMATs {} are implemented in MDAnalysis."
        "".format(", ".join(mda._READERS.keys())),
        )

    universe_group.add_argument(
        "-traj",
        dest="trajectory_format",
        type=str,
        default=None,
        help="Override automatic trajectory type detection. "
        "See trajectory for implemented formats.")

    run_group = analysis_class_parser.add_argument_group(
        title="Analysis run Parameters",
        description="Genereal parameters specific for running the analysis"
        )
    run_group.add_argument(
        "-b",
        dest="begin",
        type=str,
        default="0",
        help="frame or start time for evaluation. (default: %(default)s)"
        )

    run_group.add_argument(
        "-e",
        dest="end",
        type=str,
        default="-1",
        help="frame or end time for evaluation. (default: %(default)s)"
        )

    run_group.add_argument(
        "-dt",
        dest="dt",
        type=str,
        default="1",
        help="step or time step for evaluation. (default: %(default)s)"
        )

    run_group.add_argument(
        "-v",
        dest="verbose",
        help="Be loud and noisy",
        action="store_true"
        )

    # TODO: Should only be added if class does not have save_results 
    # function. However we only have a dict here and can not check this 
    # currently...
    save_group = analysis_class_parser.add_argument_group(
        title="Saving Parameters",
        )
    save_group.add_argument(
        "-pre",
        dest="output_prefix",
        type=str,
        default="",
        help="Additional prefix for all output files. Files will be "
             " automatically named by the used module (default: %(default)s)"
        )

    save_group.add_argument(
        "-o",
        dest="output_directory",
        type=str,
        default=".",
        help="Directory in which the output files produced will be stored."
             "(default: %(default)s)"
        )

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
        else:
            group.add_argument(
                name_par,
                dest=name,
                type=type_,
                default=default,
                help=f"{description} (default: %(default)s)"
                )
    return


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
        Customised LAMMPS `atom_style` information.

    Returns
    -------
    u : `MDAnalysis.Universe`
    """

    # MDAnalysis throws a warning if the topology does not include coordintes.
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        u = mda.Universe(topology,
                        topology_format=topology_format,
                        atom_style=atom_style)

    if coordinates is not None:
        u.load_new(coordinates,
                   format=trajectory_format)

    return u


def convert_analysis_parameters(universe,
                                analysis_callable,
                                analysis_parameters):
    """
    Convert parameters from the command line for the anlysis.

    Convert special types (i.e AtomGroups) from the command line 
    strings into the correct format.

    Parameters
    ----------
    universe : `MDAnalysis.Universe`
        Universe for which the kwargs are related to
    analysis_callable : function
        Analysis class for which the analysis is performed.
    analysis_parameters : dict
        parameters to be processed

    Returns
    -------
    converted_params : dict
        Dictionary containing the converted paramaters

    Raises
    ------
    ValueError
        If an Atomgroup does not contain any atoms
    """
    # Index 2 for paramaters
    params = parse_docs(analysis_callable)[2]
    for param_name, dictionary in params.items():
        if "AtomGroup" in dictionary['type']:
            sel = universe.select_atoms(analysis_parameters[param_name])
            if len(sel) > 0:
                analysis_parameters[param_name] = sel
            else:
                raise ValueError(
                    "AtomGroup `-{}` with selection `{}` does not "
                    "contain any atoms".format(
                        param_name,
                        analysis_parameters[param_name]
                        )
                    )
        elif "Universe" in dictionary['type']:
            analysis_parameters[param_name] = universe

    return analysis_parameters


def split_argparse_into_groups(parser, namespace):
    """
    Split the the populated namespace of argparse into groups.

    See https://stackoverflow.com/questions/31519997/is-it-possible-to-only-parse-one-argument-groups-parameters-with-argparse
    for details

    Parameters
    ----------
    parse : `argparse.ArgumentParser`
        argument parser instance
    namespace : `argparse.Namespace`
        instance storing the parameters

    Returns
    -------
    arg_grouped_dict : dict
        Dictionary containing parameters split according to their groups
    """
    arg_grouped_dict = {}
    for group in parser._action_groups:
        group_dict={a.dest:getattr(namespace, a.dest, None)
            for a in group._group_actions}
        arg_grouped_dict[group.title] = group_dict
    
    return arg_grouped_dict


def maincli(ap):
    """Execute main client interface.
    
    Returns
    -------
    ac : `MDAnalysis.analysis.base.AnalysisBase`
        AnalysisBase instance of the given ``analysis_callable`` after run.
    """
    if len(sys.argv) < 2:
        ap.error("A subcommand is required.")

    try:
        args = ap.parse_args()
        verbose = args.verbose
        analysis_callable = args.analysis_callable

        # Get the correct ArgumentParser instance from all subparsers
        # [0] selects the first subparser where our analysises live in.
        ap_sup = ap._subparsers._group_actions[0].choices[
            analysis_callable.__name__]
        arg_grouped_dict = split_argparse_into_groups(ap_sup, args)

        # Initilize Universe
        if verbose:
            print("Loading trajectory...", end="")
        universe = create_universe(**arg_grouped_dict["Universe Parameters"])
        if verbose:
            print("Done!\n")
            sys.stdout.flush()

        # Setup Analysis instance
        ac_parameters_dict = arg_grouped_dict["Mandatory Parameters"]

        # Only exsists if optional arguments are defined
        try:
            ac_parameters_dict += arg_grouped_dict["Optional Parameters"]
        except KeyError:
            pass

        ac_parameters_dict = convert_analysis_parameters(
                                universe,
                                analysis_callable,
                                ac_parameters_dict)

        ac = analysis_callable(**ac_parameters_dict)

        # Run the analysis
        with warnings.catch_warnings():
            warnings.simplefilter('always')
            startframe = convert_str_time(
                    arg_grouped_dict["Analysis run Parameters"]["begin"],
                    universe.trajectory.dt)
            stopframe = convert_str_time(
                    arg_grouped_dict["Analysis run Parameters"]["end"],
                    universe.trajectory.dt)
            step = convert_str_time(
                    arg_grouped_dict["Analysis run Parameters"]["dt"],
                    universe.trajectory.dt)

        # raises error if frame selection range is an empty selection
        if not list(range(universe.trajectory.n_frames)[slice(startframe,
                                                              stopframe,
                                                              step)]):
            raise ValueError(f"Trajectory frame range {startframe}:"
                             f"{stopframe}:{step} is not valid for"
                             f" {universe.trajectory.n_frames} frames.")

        ac.run(start=startframe,
               stop=stopframe,
               step=step,
               verbose=verbose)

        # Save the data
        output_directory = arg_grouped_dict["Saving Parameters"]["output_directory"]
        output_prefix = arg_grouped_dict["Saving Parameters"]["output_prefix"]
        output_prefix += "_" if len(output_prefix) > 0 else ""

        try:
            ac.save_results()

        except AttributeError:
            save_results(
                ac.results,
                os.path.join(
                    output_directory,
                    f"{output_prefix}{type(ac).__name__}"
                    ),
                )

    except Exception as e:
        sys.exit(Emphasise.error(f"Error: {e}"))

    return ac


def setup_clients(title, members):
    """
    Set up ArgumentParser clients.

    Parameters
    ----------
    title : str
        title of the parser
    members : list
        list containing Analysis classes for setting up the parser

    Returns
    -------
    argparse.ArgumentParser instance
    """
    ap = argparse.ArgumentParser()
    ap.add_argument('--version',
                    action='version',
                    version="mdacli {}".format(__version__))

    cli_parser = ap.add_subparsers(title=title)

    analysis_interfaces = {
        member.__name__: parse_callable_signature(member)
        for member in members
        }

    # adds each Analysis class/function as a CLI under 'cli_parser'
    # to be writen
    for interface_name, parameters in analysis_interfaces.items():
        create_CLI(cli_parser, interface_name, parameters)

    return ap


def main():
    """Execute main CLI entry point."""
    members = find_AnalysisBase_members_ignore_warnings(_relevant_modules)

    if members is None:
        sys.exit("No analysis modules found.")

    title = "MDAnalysis Analysis CLI"
    ap = setup_clients(title=title, members=members)
    maincli(ap)
