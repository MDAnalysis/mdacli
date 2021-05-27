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
from MDAnalysisTests.datafiles import TPR, XTC
from MDAnalysisTests.topology.test_lammpsdata import LAMMPS_NORESID
from MDAnalysisTests.core.test_universe import CHOL_GRO

from mdacli.cli import _relevant_modules, setup_clients
from mdacli.cli import create_universe

from mdacli.libcli import find_AnalysisBase_members_ignore_warnings


_relev_mod = list(_relevant_modules)


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
     ('-f', "coordinates", ["foo", "bar"]),
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
    members = find_AnalysisBase_members_ignore_warnings(_relev_mod)

    with patch.object(sys, 'argv', testargs):
        args = setup_clients(title='title', members=members).parse_args()
        t = type(val)
        assert t(getattr(args, dest)) == val


class Test_create_Universe:
    def test_create_universe(self):
        u = create_universe(TPR, XTC)
        assert len(u.atoms) == 47681
        assert hasattr(u, "trajectory")

    def test_topology_format(self):
        u = create_universe(StringIO(CHOL_GRO), topology_format="GRO")
        assert len(u.atoms) == 8

    def test_create_universe_no_traj(self):
        u = create_universe(TPR)
        assert not hasattr(u, "trajectory")

    def test_trajectory_format(self):
        u = create_universe(topology=StringIO(CHOL_GRO),
                            coordinates=StringIO(CHOL_GRO),
                            topology_format="GRO",
                            trajectory_format="GRO")
        assert len(u.atoms) == 8

    def test_create_universe_atom_style(self):
        u = create_universe(topology=StringIO(LAMMPS_NORESID),
                            topology_format="data",
                            atom_style="id type x y z")
        assert u.atoms[0].mass ==  28.0


# class Test_run_analsis(object):
#     """Test class for analyze_data."""

#     @pytest.fixture()
#     def kwargs(self):
#         """Keyword arguments for run."""
#         kwargs = {}
#         kwargs["topology"] = TPR
#         kwargs["topology_format"] = None
#         kwargs["trajectories"] = XTC
#         kwargs["trajectory_format"] = None
#         kwargs["atom_style"] = None
#         kwargs["begin"] = "0"
#         kwargs["end"] = "1"
#         kwargs["dt"] = "1"
#         kwargs["output_directory"] = "."
#         kwargs["output_prefix"] = ""
#         kwargs["verbose"] = False
#         kwargs["g1"] = "resid 500"  # water molecules
#         kwargs["g2"] = "resid 501"  # water molecules

#         kwargs["func"] = None
#         return kwargs

#     def test_default(self, kwargs, tmpdir):
#         """Simple with default arguments."""
#         with tmpdir.as_cwd():
#             run_analsis(analysis_callable=InterRDF, **kwargs)

#     def test_topology_format(self, kwargs, tmpdir):
#         """Simple test run."""
#         kwargs["topology_format"] = "TPR"
#         with tmpdir.as_cwd():
#             run_analsis(analysis_callable=InterRDF, **kwargs)

#     def test_atom_style(self, kwargs, tmpdir):
#         """Test the LAMMPS data format."""
#         kwargs["topology"] = StringIO(LAMMPS_NORESID)
#         kwargs["topology_format"] = "data"
#         kwargs["atom_style"] = "id type x y z"
#         kwargs["trajectories"] = None
#         kwargs["g1"] = "all"
#         kwargs["g2"] = "all"

#         with tmpdir.as_cwd():
#             ac = run_analsis(analysis_callable=InterRDF, **kwargs)

#         assert len(ac.u.atoms) == 1
#         assert_equal(ac.u.atoms[0].mass, 28.0)

#     def test_trajectory_format(self, kwargs, tmpdir):
#         """Test for trajectory format."""
#         kwargs["trajectory_format"] = "XTC"
#         with tmpdir.as_cwd():
#             ac = run_analsis(analysis_callable=InterRDF, **kwargs)
#             assert ac._trajectory.format == "XTC"

#     def test_verbose(self, capsys, kwargs, tmpdir):
#         """Test for being verbose."""
#         kwargs["verbose"] = True
#         with tmpdir.as_cwd():
#             run_analsis(analysis_callable=InterRDF, **kwargs)
#         captured = capsys.readouterr()
#         assert "Loading trajectory..." in captured.out

#     def test_traj_clices(self, kwargs, tmpdir):
#         """Test for trajectory slicing."""
#         kwargs["begin"] = "1"
#         kwargs["end"] = "10"
#         kwargs["dt"] = "2"
#         with tmpdir.as_cwd():
#             ac = run_analsis(analysis_callable=InterRDF, **kwargs)
#         assert ac.start == 1
#         assert ac.stop == 10
#         assert ac.step == 2
