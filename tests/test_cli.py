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
from io import StringIO
from unittest.mock import patch

import pytest
from MDAnalysis.analysis.rdf import InterRDF
from MDAnalysisTests.datafiles import DCD, PSF
from MDAnalysisTests.topology.test_lammpsdata import LAMMPS_NORESID
from numpy.testing import assert_equal

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
     ('-dt', "dt", 42)))
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


class Test_run_analsis(object):
    """Test class for analyze_data."""

    @pytest.fixture()
    def kwargs(self):
        """Keyword arguments for run."""
        kwargs = {}
        kwargs["topology"] = PSF
        kwargs["topology_format"] = None
        kwargs["trajectories"] = DCD
        kwargs["trajectory_format"] = None
        kwargs["atom_style"] = None
        kwargs["begin"] = "0"
        kwargs["end"] = "1"
        kwargs["dt"] = "1"
        kwargs["verbose"] = False
        kwargs["g1"] = "all"
        kwargs["g2"] = "all"

        kwargs["func"] = None
        return kwargs

    def test_default(self, kwargs, tmpdir):
        """Simple with default arguments."""
        with tmpdir.as_cwd():
            run_analsis(analysis_callable=InterRDF, **kwargs)

    def test_topology_format(self, kwargs, tmpdir):
        """Simple test run."""
        kwargs["topology_format"] = "PSF"
        with tmpdir.as_cwd():
            run_analsis(analysis_callable=InterRDF, **kwargs)

    def test_atom_style(self, kwargs, tmpdir):
        """Test the LAMMPS data format."""
        kwargs["topology"] = StringIO(LAMMPS_NORESID)
        kwargs["topology_format"] = "data"
        kwargs["atom_style"] = "id type x y z"
        kwargs["trajectories"] = None

        with tmpdir.as_cwd():
            ac = run_analsis(analysis_callable=InterRDF, **kwargs)

        assert len(ac.u.atoms) == 1
        assert_equal(ac.u.atoms[0].mass, 28.0)

    def test_trajectory_format(self, kwargs, tmpdir):
        """Test for trajectory format."""
        kwargs["trajectory_format"] = "DCD"
        with tmpdir.as_cwd():
            ac = run_analsis(analysis_callable=InterRDF, **kwargs)
            assert ac._trajectory.format == "DCD"

    def test_verbose(self, capsys, kwargs, tmpdir):
        """Test for being verbose."""
        kwargs["verbose"] = True
        with tmpdir.as_cwd():
            run_analsis(analysis_callable=InterRDF, **kwargs)
        captured = capsys.readouterr()
        assert "Loading trajectory..." in captured.out

    def test_traj_clices(self, kwargs, tmpdir):
        """Test for trajectory slicing."""
        kwargs["begin"] = "1"
        kwargs["end"] = "10"
        kwargs["dt"] = "2"
        with tmpdir.as_cwd():
            ac = run_analsis(analysis_callable=InterRDF, **kwargs)
        assert ac.start == 1
        assert ac.stop == 10
        assert ac.step == 2
