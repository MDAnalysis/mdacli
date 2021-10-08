#!/usr/bin/env python3
# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
#
# Copyright (c) 2021 Authors and contributors
#
# Released under the GNU Public Licence, v2 or any higher version
# SPDX-License-Identifier: GPL-2.0-or-later

"""Functionalities that support the creation of the command interfaces."""
import argparse
import importlib
import inspect
import json
import warnings

import MDAnalysis as mda
from MDAnalysis.analysis.base import AnalysisBase


# UNIVERSES = defaultdict(dict)


def store_universe(name):
    class StoreUniverse(argprase.Action):
        def __call__(self, parser, namespace, value, option_string=None):
            UNIVERSES[name][self.dest] = value
    return StoreUniverse


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


def find_AnalysisBase_members(modules):
    """Find Analysis Base members in modules."""
    members = find_classes_in_modules(
        AnalysisBase, *[f'MDAnalysis.analysis.{m}' for m in modules])
    return members


def find_AnalysisBase_members_ignore_warnings(modules):
    """Find Analysis Base members in modules and ignoring all warnings."""
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        members = find_AnalysisBase_members(modules)
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


def add_analysis_run_parameters(analysis_class_parser):
    run_group = analysis_class_parser.add_argument_group(
        title="Analysis Run Parameters",
        description="Genereal parameters specific for running the analysis"
        )
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
        action="store_true"
        )

def add_output_group(analysis_class_parser):
    output_group = analysis_class_parser.add_argument_group(
        title="Output Parameters",
        )
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
    """."""
    name = f'-{name}' if name else ''
    parser.add_argument(
        f"-s{name}",
        dest="topology",
        type=str,
        default="topol.tpr",
        action=store_universe(name),
        help="The topolgy file. "
        "The FORMATs {} are implemented in MDAnalysis."
        "".format(", ".join(mda._PARSERS.keys())),
        )

    parser.add_argument(
        f"-top{name}",
        dest="topology_format",
        type=str,
        action=store_universe(name),
        default=None,
        help="Override automatic topology type detection. "
        "See topology for implemented formats.")

    parser.add_argument(
        f"-atom_style{name}",
        dest="atom_style",
        type=str,
        action=store_universe(name),
        default=None,
        help="Manually set the atom_style information"
        "(currently only LAMMPS parser). E.g. atom_style='id type x y z'.")

    parser.add_argument(
        f"-f{name}",
        dest="coordinates",
        type=str,
        default=None,
        action=store_universe(name),
        nargs="+",
        help="A single or multiple coordinate files. "
        "The FORMATs {} are implemented in MDAnalysis."
        "".format(", ".join(mda._READERS.keys())),
        )

    parser.add_argument(
        f"-traj{name}",
        dest="trajectory_format",
        type=str,
        action=store_universe(name),
        default=None,
        help="Override automatic trajectory type detection. "
        "See trajectory for implemented formats.")

    return


def add_cli_atom_group(parser, name=''):
    """."""
    add_cli_universe(parser, name=name)

    parser.add_argument(
        '-selection',
        dest="selection",
        type=str,
        default="all",
        help="Atom selection to create MDAnalysis AtomGroup.",
        )

    return


def add_cli_single_atom_group():
    pass


def add_cli_single_universe():
    pass


def add_cli_several_universes():
    pass


def add_cli_several_atomgroups():
    pass
