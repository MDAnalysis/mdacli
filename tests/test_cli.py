#!/usr/bin/env python3
# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
#
# Copyright (c) 2021 Authors and contributors
#
# Released under the GNU Public Licence, v2 or any higher version
# SPDX-License-Identifier: GPL-2.0-or-later
"""Test mdacli cli."""
import argparse
import os
import subprocess
import sys
from io import StringIO
from unittest.mock import patch

import pytest
from MDAnalysis.analysis.gnm import GNMAnalysis
from MDAnalysis.analysis.rdf import InterRDF
from MDAnalysis.analysis.rms import RMSF
from MDAnalysis.core.universe import Universe
from MDAnalysisTests.core.test_universe import CHOL_GRO
from MDAnalysisTests.datafiles import TPR, XTC
from MDAnalysisTests.topology.test_lammpsdata import LAMMPS_NORESID

from mdacli.cli import (
    _relevant_modules,
    convert_analysis_parameters,
    create_universe,
    run_analsis,
    setup_clients,
    )
from mdacli.libcli import find_AnalysisBase_members_ignore_warnings


_relev_mod = list(_relevant_modules)


def test_required_args():
    """Test that there is a module given."""
    with pytest.raises(subprocess.CalledProcessError):
        subprocess.check_call(['mda'])


def test_wrong_module():
    """Test for a non existent module."""
    with pytest.raises(subprocess.CalledProcessError):
        subprocess.check_call(['mda', 'foo'])


@pytest.mark.parametrize('args', ("version", "help"))
def test_extra_options(args):
    """Test for a ab extra option."""
    subprocess.check_call(['mda', '--' + args])


@pytest.mark.parametrize('args', ("RMSF", "rmsf"))
def test_case_insensitive(args):
    """Test for beeing case insensitive."""
    subprocess.check_call(['mda', args, "-h"])


@pytest.mark.parametrize(
    'opt, dest, val',
    (('-s', "topology", "foo"),
     ('-top', "topology_format", "foo"),
     ('-f', "coordinates", ["foo", "bar"]),
     ('-traj', "trajectory_format", "foo"),
     ('-atom_style', "atom_style", "foo"),
     ('-b', "start", 42),
     ('-e', "stop", 42),
     ('-dt', "step", 42)))
def test_setup_clients(opt, dest, val):
    """Test all additional arguments."""
    testargs = ["mdacli", "rmsf", opt]
    if type(val) == list:
        for i in val:
            testargs.append(str(i))
    else:
        testargs.append(str(val))
    members = find_AnalysisBase_members_ignore_warnings(_relev_mod)

    with patch.object(sys, 'argv', testargs):
        ap = argparse.ArgumentParser()
        setup_clients(ap, title='title', members=members)
        args = ap.parse_args()
        t = type(val)
        assert t(getattr(args, dest)) == val


class Test_create_Universe:
    """Test initilizing mda universes."""

    def test_create_universe(self):
        """Test creating a simple Universe."""
        u = create_universe(TPR, XTC)
        assert len(u.atoms) == 47681
        assert hasattr(u, "trajectory")

    def test_topology_format(self):
        """Test non default topology format."""
        u = create_universe(StringIO(CHOL_GRO), topology_format="GRO")
        assert len(u.atoms) == 8

    def test_create_universe_no_traj(self):
        """Test non trajectpry parsed."""
        u = create_universe(TPR)
        assert not hasattr(u, "trajectory")

    def test_trajectory_format(self):
        """Test non default trajectpry format."""
        u = create_universe(topology=StringIO(CHOL_GRO),
                            coordinates=StringIO(CHOL_GRO),
                            topology_format="GRO",
                            trajectory_format="GRO")
        assert len(u.atoms) == 8

    def test_create_universe_atom_style(self):
        """Test custom atom style."""
        u = create_universe(topology=StringIO(LAMMPS_NORESID),
                            topology_format="data",
                            atom_style="id type x y z")
        assert u.atoms[0].mass == 28.0


class Test_convert_analysis_parameters:
    """Test class for converting analysis parameters."""

    @pytest.fixture()
    def universe(self):
        """Univserse fixture."""
        return Universe(TPR, XTC)

    def test_Atomgroup(self, universe):
        """Test AtomGroup conversion."""
        analysis_parameters = {"atomgroup": "all"}
        convert_analysis_parameters(universe,
                                    RMSF,
                                    analysis_parameters)
        assert analysis_parameters["atomgroup"] == universe.atoms

    def test_Universe(self, universe):
        """Test Universe conversion."""
        analysis_parameters = {"universe": None}
        convert_analysis_parameters(universe,
                                    GNMAnalysis,
                                    analysis_parameters)
        assert analysis_parameters["universe"] == universe

    def test_zero_atoms(self, universe):
        """Test error if zero atoms are present after conversion."""
        analysis_parameters = {"atomgroup": "name foo"}
        with pytest.raises(ValueError, match="AtomGroup `-atomgroup` with "):
            convert_analysis_parameters(universe,
                                        RMSF,
                                        analysis_parameters)

    def test_only_set_if_key_exists(self, universe):
        """
        Test if keys exist.

        The docparser does not distinguis between mandatory and positional
        arguments. So we only should set the value if it is inside the dict.
        """
        analysis_parameters = {}
        convert_analysis_parameters(universe,
                                    RMSF,
                                    analysis_parameters)
        assert not analysis_parameters


class Test_run_analsis:
    """Test class for analyze_data."""

    @pytest.fixture()
    def universe_parameters(self):
        """Universe fixture."""
        kwargs = {}
        kwargs["topology"] = TPR
        kwargs["topology_format"] = None
        kwargs["coordinates"] = XTC
        kwargs["trajectory_format"] = None
        kwargs["atom_style"] = None
        return kwargs

    @pytest.fixture()
    def mandatory_parameters(self):
        """Mandatory parameter fixture."""
        kwargs = {}
        kwargs["g1"] = "resid 500"  # water molecules
        kwargs["g2"] = "resid 501"  # water molecules
        return kwargs

    def test_default(self,
                     universe_parameters,
                     mandatory_parameters,
                     tmpdir):
        """Test with default arguments."""
        with tmpdir.as_cwd():
            run_analsis(analysis_callable=InterRDF,
                        universe_parameters=universe_parameters,
                        mandatory_analysis_parameters=mandatory_parameters)

    def test_optional_paramaters(self,
                                 universe_parameters,
                                 mandatory_parameters,
                                 tmpdir):
        """Test with optional parameters given."""
        opt_params = {"nbins": 100}
        with tmpdir.as_cwd():
            a = run_analsis(analysis_callable=InterRDF,
                            universe_parameters=universe_parameters,
                            mandatory_analysis_parameters=mandatory_parameters,
                            optional_analysis_parameters=opt_params)

        assert a.rdf_settings["bins"] == opt_params["nbins"]

    def test_verbose(self,
                     universe_parameters,
                     mandatory_parameters,
                     tmpdir,
                     capsys):
        """Test for being verbose."""
        run_parameters = {"verbose": True}
        with tmpdir.as_cwd():
            run_analsis(analysis_callable=InterRDF,
                        universe_parameters=universe_parameters,
                        mandatory_analysis_parameters=mandatory_parameters,
                        run_parameters=run_parameters)
        captured = capsys.readouterr()
        assert "Loading trajectory..." in captured.out

    def test_custom_output(self,
                           universe_parameters,
                           mandatory_parameters,
                           tmpdir):
        """Test for custom output."""
        output_parameters = {"output_directory": "foo",
                             "output_prefix": "bar"}

        with tmpdir.as_cwd():
            os.mkdir("foo")
            run_analsis(analysis_callable=InterRDF,
                        universe_parameters=universe_parameters,
                        mandatory_analysis_parameters=mandatory_parameters,
                        output_parameters=output_parameters)

            os.path.isfile(os.path.join("foo",
                                        "bar_InterRDF_count_bins_rdf.csv"))

    def test_output_directory(self,
                              universe_parameters,
                              mandatory_parameters,
                              tmpdir):
        """Test for custom output directory."""
        output_parameters = {"output_directory": "foo"}
        with tmpdir.as_cwd():
            os.mkdir("foo")
            run_analsis(analysis_callable=InterRDF,
                        universe_parameters=universe_parameters,
                        mandatory_analysis_parameters=mandatory_parameters,
                        output_parameters=output_parameters)

            os.path.isfile(os.path.join("foo", "InterRDF_count_bins_rdf.csv"))

    def test_output_prefix(self,
                           universe_parameters,
                           mandatory_parameters,
                           tmpdir):
        """Test for custom prefix."""
        output_parameters = {"output_prefix": "foo"}
        with tmpdir.as_cwd():
            os.mkdir("foo")
            run_analsis(analysis_callable=InterRDF,
                        universe_parameters=universe_parameters,
                        mandatory_analysis_parameters=mandatory_parameters,
                        output_parameters=output_parameters)

            os.path.isfile("foo_InterRDF_count_bins_rdf.csv")
