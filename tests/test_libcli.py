#!/usr/bin/env python3
# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
#
# Copyright (c) 2020 Authors and contributors
#
# Released under the GNU Public Licence, v2 or any higher version
# SPDX-License-Identifier: GPL-2.0-or-later
"""Test libcli."""
import argparse
import os
from json.decoder import JSONDecodeError

import pytest
from MDAnalysis.analysis.base import AnalysisBase
from MDAnalysis.analysis.helix_analysis import HELANAL
from MDAnalysis.analysis.lineardensity import LinearDensity
from MDAnalysis.analysis.rdf import InterRDF, InterRDF_s

from mdacli import libcli

from . import example_json


@pytest.mark.parametrize(
    'cmd,expected',
    [
        ('-d {"key1":1}', {"key1": 1}),
        ]
    )
def test_KwargsDict(cmd, expected):
    """Test dict reading action."""
    ap = argparse.ArgumentParser()
    ap.add_argument(
        '-d',
        action=libcli.KwargsDict,
        default=None,
        )
    args = ap.parse_args(cmd.split())
    assert args.d == expected


@pytest.mark.parametrize(
    'cmd,expected',
    [
        (f"-d {os.fspath(example_json)}", {"key1": 1}),
        ]
    )
def test_KwargsDict_from_file(cmd, expected):
    """Test dict reading action."""
    ap = argparse.ArgumentParser()
    ap.add_argument(
        '-d',
        action=libcli.KwargsDict,
        default=None,
        )
    args = ap.parse_args(cmd.split())
    assert args.d == expected


def test_KwargsDict_error():
    """Test error."""
    ap = argparse.ArgumentParser()
    ap.add_argument(
        '-d',
        action=libcli.KwargsDict,
        default=None,
        )
    with pytest.raises(JSONDecodeError):
        ap.parse_args("-d fail".split())


def test_find_AnalysisBase_members():
    """Test several input modules."""
    names = ["helix_analysis", "lineardensity"]
    members = libcli.find_AnalysisBase_members(names)
    assert members[0] is HELANAL
    assert members[1] is LinearDensity


def test_find_AnalysisBase_members_single():
    """Test one input module."""
    members = libcli.find_AnalysisBase_members(['rdf'])

    assert members[0] is InterRDF
    assert members[1] is InterRDF_s


def test_find_AnalysisBase_members_None():
    """Test that no module is found."""
    members = libcli.find_classes_in_modules(AnalysisBase, "MDAnalysis")

    assert members is None


def split_argparse_into_groups():
    """Test splitting argparse Namespace into several dicts."""
    parser = argparse.ArgumentParser()

    group1 = parser.add_argument_group('group1')
    group1.add_argument('--test1', help="test1")

    group2 = parser.add_argument_group('group2')
    group2.add_argument('--test2', help="test2")

    args = parser.parse_args('--test1 one --test2 two'.split())

    arg_grouped_dict = split_argparse_into_groups(parser, args)

    assert arg_grouped_dict["group1"]["test1"] == "one"
    assert arg_grouped_dict["group2"]["test2"] == "two"
