#!/usr/bin/env python3
# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
#
# Copyright (c) 2021 Authors and contributors
#
# Released under the GNU Public Licence, v2 or any higher version
# SPDX-License-Identifier: GPL-2.0-or-later

"""Functionalities that support the creation of the command lines interface."""
import argparse
import importlib
import inspect
import json
import logging
import os
import warnings

import MDAnalysis as mda
from MDAnalysis.analysis.base import AnalysisBase
from MDAnalysis.transformations.boxdimensions import set_dimensions

from .colors import Emphasise
from .save import save_results
from .utils import convert_str_time, parse_callable_signature, parse_docs


logger = logging.getLogger(__name__)

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


class KwargsDict(argparse.Action):
    """
    Convert input string to a dictionary.

    If string points to a ".json" file, reads the file.
    Else, attempts to convert string to dictionary using json.loads.
    """

    def __call__(self, parser, namespace, value, option_string=None):
        """Call me."""
        if value.startswith("{") and value.endswith("}"):
            try:
                jdict = json.loads(value)

            except json.decoder.JSONDecodeError as err:
                raise json.decoder.JSONDecodeError(
                    "An error ocurred when reading "
                    f"{self.dest!r} argument.",
                    err.doc,
                    err.pos,
                    ) from None
        else:
            with open(value, 'r') as fin:
                jdict = json.load(fin)

        setattr(namespace, self.dest, jdict)


def find_classes_in_modules(klass, *module_names):
    """
    Find classes in modules.

    A series of names can be given as arguments.

    Parameters
    ----------
    klass : class type
        The class type to search for.

    module_names : str
        module to import import in absolute or relative terms
        (e.g. either pkg.mod or ..mod).

    Returns
    -------
    list of found class objects.
    If no classes are found, return None.
    """
    members = []
    for name in module_names:
        module = importlib.import_module(name)
        for _, member in inspect.getmembers(module):
            if inspect.isclass(member) and issubclass(member, klass) \
               and member is not klass:
                members.append(member)

    return members or None


def find_AnalysisBase_members(modules,
                              ignore_warnings=False):
    """Find Analysis Base members in modules.

    Parameters
    ----------
    modules : list
        list of modules for which members should be searched for
    ignore_warnings : bool
        Flag to ignore warnings
    """
    with warnings.catch_warnings():
        if not ignore_warnings:
            warnings.simplefilter('ignore')
        members = find_classes_in_modules(
            AnalysisBase, *[f'MDAnalysis.analysis.{m}' for m in modules])
    return members


def split_argparse_into_groups(parser, namespace):
    """
    Split the populated namespace of argparse into groups.

    https://stackoverflow.com/questions/31519997/is-it-possible-to-
    only-parse-one-argument-groups-parameters-with-argparse

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
        group_dict = {a.dest: getattr(namespace, a.dest, None)
                      for a in group._group_actions}
        arg_grouped_dict[group.title] = group_dict

    return arg_grouped_dict


def add_run_group(analysis_class_parser):
    """Add run group parameters to an given argparse.ArgumentParser instance.

    The run group adds the parameters `start`, `stop`, `step`, `verbose` to the
    parser.

    Parameters
    ----------
    analysis_class_parser : argparse.ArgumentParser
        The ArgumentsParser instance to which the run grorup is added
    """
    run_group = analysis_class_parser.add_argument_group(
        title="Analysis Run Parameters",
        description="Genereal parameters specific for running the analysis")

    run_group.add_argument(
        "-b",
        dest="start",
        type=str,
        default="0",
        help="frame or start time for evaluation. (default: %(default)s)"
        )

    run_group.add_argument(
        "-e",
        dest="stop",
        type=str,
        default="-1",
        help="frame or end time for evaluation. (default: %(default)s)"
        )

    run_group.add_argument(
        "-dt",
        dest="step",
        type=str,
        default="1",
        help="step or time step for evaluation. (default: %(default)s)"
        )

    run_group.add_argument(
        "-v",
        dest="verbose",
        help="Be loud and noisy",
        action="store_true",
        )


def add_output_group(analysis_class_parser):
    """Add output group parameters to argparse.ArgumentParser instance.

    The run group adds the parameters `output_prefix` and
    `output_directory` to the parser.

    Parameters
    ----------
    analysis_class_parser : argparse.ArgumentParser
        The ArgumentsParser instance to which the run grorup is added
    """
    output_group = analysis_class_parser.add_argument_group(
        title="Output Parameters",
        description="Genereal parameters specific for the result output.")

    output_group.add_argument(
        "-pre",
        dest="output_prefix",
        type=str,
        default="",
        help="Additional prefix for all output files. Files will be "
             " automatically named by the used module (default: %(default)s)"
        )

    output_group.add_argument(
        "-o",
        dest="output_directory",
        type=str,
        default=".",
        help="Directory in which the output files produced will be stored."
             "(default: %(default)s)"
        )


def add_cli_universe(parser, name=''):
    """Add universe parameters to an given argparse.ArgumentParser.

    instance. The parameters `topology`, `topology_format`, `atom_style`,
    `coordinates` and `trajectory_format` are added to the parse.

    Parameters
    ----------
    analysis_class_parser : argparse.ArgumentParser
        The ArgumentsParser instance to which the run grorup is added
    name : str
        suffix for the argument names
    """
    name = f'_{name}' if name else ''

    parser.add_argument(
        f"-s{name}",
        dest=f"topology{name}",
        type=str,
        default="topol.tpr",
        help="The topolgy file. "
        "The FORMATs {} are implemented in MDAnalysis."
        "".format(", ".join(mda._PARSERS.keys())),
        )

    parser.add_argument(
        f"-top{name}",
        dest=f"topology_format{name}",
        type=str,
        default=None,
        help="Override automatic topology type detection. "
        "See topology for implemented formats.")

    parser.add_argument(
        f"-atom_style{name}",
        dest=f"atom_style{name}",
        type=str,
        default=None,
        help="Manually set the atom_style information"
        "(currently only LAMMPS parser). E.g. atom_style='id type x y z'.")

    parser.add_argument(
        f"-f{name}",
        dest=f"coordinates{name}",
        type=str,
        default=None,
        nargs="+",
        help="A single or multiple coordinate files. "
        "The FORMATs {} are implemented in MDAnalysis."
        "".format(", ".join(mda._READERS.keys())),
        )

    parser.add_argument(
        f"-traj{name}",
        dest=f"trajectory_format{name}",
        type=str,
        default=None,
        help="Override automatic trajectory type detection. "
        "See trajectory for implemented formats.")

    parser.add_argument(
        f"-dimensions{name}",
        dest=f"dimensions{name}",
        type=float,
        default=None,
        nargs="+",
        help="Manually set/overwrite the simulation box dimensions to a "
        "vector containing unit cell dimensions [a, b, c, α, β, γ], "
        "lengths a, b, c are in Å, and angles α, β, γ are in degrees. "
        "Providing only three parameters will assume a rectengular simulation "
        "box (α = β = γ = 90°).")


def create_CLI(sub_parser, interface_name, parameters):
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
    sub_parser : argparse.sub_parser
        A sub parser where the new parser will be added.

    interface_name : str
        Name of the interface.

    parameters : dict
        Parameters needed to fill the argparse requirements for the
        CLI interface.

    Returns
    -------
    None
    """
    # creates the subparser
    analysis_class_parser = sub_parser.add_parser(
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


def create_universe(topology,
                    coordinates=None,
                    topology_format=None,
                    trajectory_format=None,
                    atom_style=None,
                    dimensions=None):
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
    dimensions : iterable of floats
        vector that contains unit cell lengths and probable angles.
        Expected shapes are eithere (6, 0) or (1, 6) or for
        shapes of (3, 0) or (1, 3) all angles are set to 90 degrees.

    Raises
    ------
    IndexError
        If the dimesions of the `dimensions` argument are not 3 or 6.

    Returns
    -------
    `MDAnalysis.Universe`
    """
    universe = mda.Universe(topology,
                            topology_format=topology_format,
                            atom_style=atom_style)

    if coordinates is not None:
        universe.load_new(coordinates, format=trajectory_format)

    if dimensions is not None:
        if len(dimensions) == 3:
            dimensions = [*dimensions, 90, 90, 90]
        elif len(dimensions) != 6:
            raise IndexError(
                "The dimensions must contain at least 3 entries for "
                "the box length and possibly 3 more entries for the angles.")

        universe.trajectory.add_transformations(set_dimensions(dimensions))

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
    analysis_parameters_keys = [p.split("_")[-1] for p
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
                # Create universe parameter dictionary from signature
                sig = inspect.signature(create_universe)
                universe_parameters = dict(sig.parameters)

                for k in universe_parameters.keys():
                    universe_parameters[k] = analysis_parameters.pop(
                        f"{k}_{param_name}")

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
    for member_name, parameters in analysis_interfaces.items():
        create_CLI(sub_parser=cli_subparser,
                   interface_name=member_name.lower(),
                   parameters=parameters)


def init_base_argparse(name, version, description):
    """Create a basic `ArgumentParser`.

    The parser has options for printing the version, running in debug mode
    and with a logfile. Note that the funtion only adds the options to
    the parser but not the logic for actually running in debug mode nor
    how to store the log file.

    Parameters
    ----------

    name : str
        Name of the cli program

    version : str
        Version of the cli program

    description : str
        Description of the cli program

    Returns
    -------
    `ArgumentParser`
    """
    ap = argparse.ArgumentParser(
        description=description,
        formatter_class=argparse.RawDescriptionHelpFormatter)

    ap.add_argument(
        '--version',
        action='version',
        version=f"{name} {version}",
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
