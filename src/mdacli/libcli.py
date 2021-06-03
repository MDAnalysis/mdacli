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
from pathlib import Path

from MDAnalysis.analysis.base import AnalysisBase


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
    Split the the populated namespace of argparse into groups.

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
