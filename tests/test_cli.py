#!/usr/bin/env python3
# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
#
# Copyright (c) 2021 Authors and contributors
#
# Released under the GNU Public Licence, v2 or any higher version
# SPDX-License-Identifier: GPL-2.0-or-later
"""Test mdacli cli."""
import subprocess
import sys
from unittest.mock import patch

import pytest
from MDAnalysis.analysis.rdf import InterRDF
from MDAnalysisTests.datafiles import DCD, PSF

from mdacli.cli import run_analsis, setup_clients


def test_required_args():
    """Test that there is a module given."""
    with pytest.raises(subprocess.CalledProcessError):
        subprocess.check_call(['mdacli'])


def test_wrong_module():
    """Test for a non existent module."""
    with pytest.raises(subprocess.CalledProcessError):
        subprocess.check_call(['mdacli', 'foo'])


@pytest.mark.parametrize('args', ("version", "help"))
def test_extra_options(args):
    """Test for a ab extra option."""
    subprocess.check_call(['mdacli', '--' + args])


@pytest.mark.parametrize(
    'opt, dest, val',
    (('-s', "topology", "foo"),
     ('-top', "topology_format", "foo"),
     ('-f', "trajectories", ["foo", "bar"]),
     ('-traj', "trajectory_format", "foo"),
     ('-atom_style', "atom_style", "foo"),
     ('-b', "begin", 42),
     ('-e', "end", 42),
     ('-dt', "dt", 42),
     ('-box', "box", [42, 42, 42])))
def test_setup_clients(opt, dest, val):
    """Test all additional arguments."""
    testargs = ["mdacli", "RMSF", opt]
    if type(val) == list:
        for i in val:
            testargs.append(str(i))
    else:
        testargs.append(str(val))
    with patch.object(sys, 'argv', testargs):
        args = setup_clients().parse_args()
        t = type(val)
        assert t(getattr(args, dest)) == val


class Test_analyze_data(object):
    """Test class for analyze_data."""

    @pytest.fixture()
    def kwargs(self):
        """Keyword arguments for run."""
        kwargs = {}
        kwargs["topology"] = PSF
        kwargs["trajectories"] = DCD
        kwargs["begin"] = "0"
        kwargs["begin"] = "0"
        kwargs["end"] = "1"
        kwargs["dt"] = "1"
        kwargs["verbose"] = False
        kwargs["g1"] = "all"
        kwargs["g2"] = "all"
        kwargs["box"] = None
        kwargs["atom_style"] = None
        kwargs["trajectory_format"] = None
        kwargs["topology_format"] = None
        kwargs["func"] = None
        return kwargs

    def test_run(self, kwargs, tmpdir):
        """Simple test run."""
        with tmpdir.as_cwd():
            run_analsis(analysis_callable=InterRDF, **kwargs)
