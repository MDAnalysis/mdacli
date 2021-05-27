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
        try:
            if Path(value).exists():
                # read dict from JSON file
                with open(value, 'r') as fin:
                    jdict = json.load(fin)

            else:
                jdict = json.loads(value)

        except json.decoder.JSONDecodeError as err:
            raise json.decoder.JSONDecodeError(
                "An error ocurred when reading "
                f"{self.dest!r} argument.",
                err.doc,
                err.pos,
                ) from None

        setattr(namespace, self.dest, jdict)


def find_AnalysisBase_members(*module_names):
    """
    Check for AnalysiBase members in the module given by module_names.

    A series of names can be given as arguments.

    Parameters
    ----------
    module_names : str
        module to import import in absolute or relative terms
        (e.g. either pkg.mod or ..mod).

    Returns
    -------
    list of analysis classes
    """
    members = []
    for name in module_names:
        module = importlib.import_module(name)
        for _, member in inspect.getmembers(module):
            if inspect.isclass(member) and issubclass(member, AnalysisBase) \
               and member is not AnalysisBase:
                members.append(member)

    return None if not members else members
