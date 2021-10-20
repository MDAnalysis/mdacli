#!/usr/bin/env python3
# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
#
# Copyright (c) 2021 Authors and contributors
#
# Released under the GNU Public Licence, v2 or any higher version
# SPDX-License-Identifier: GPL-2.0-or-later
"""Test mdacli cli."""
import argparse
import logging
import os
import subprocess
import sys
from io import StringIO
from unittest.mock import patch

import pytest
from MDAnalysis.analysis.msd import EinsteinMSD
from MDAnalysis.analysis.rdf import InterRDF
from MDAnalysis.analysis.rms import RMSF
from MDAnalysis.core.universe import Universe
from MDAnalysisTests.core.test_universe import CHOL_GRO
from MDAnalysisTests.datafiles import TPR, XTC
from MDAnalysisTests.topology.test_lammpsdata import LAMMPS_NORESID

from mdacli.cli import (
    _exit_if_a_is_b,
    _relevant_modules,
    convert_analysis_parameters,
    create_universe,
    init_base_argparse,
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


@pytest.mark.parametrize('args', ("version", "debug", "help"))
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


def test__exit_if_a_is_b():
    """Test for a SystemExit using pytest.raises."""
    msg = "foo"
    with pytest.raises(SystemExit, match=msg) as error:
        _exit_if_a_is_b(1, 1, msg=msg)
    assert error.type == SystemExit


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


class Test_init_base_argparse():
    """Test for basic argument parser."""

    @pytest.fixture()
    def ap(self):
        """Return the basic parser."""
        return init_base_argparse(name="foo",
                                  version="0.0.0",
                                  description="bar")

    def test_version(self, ap, capsys):
        """Test version option."""
        with pytest.raises(SystemExit) as error:
            ap.parse_args(["--version"])

        assert error.type == SystemExit

        captured = capsys.readouterr()
        assert "foo 0.0.0\n" == captured.out

    def test_description(self, ap, capsys):
        """Test version option."""
        with pytest.raises(SystemExit) as error:
            ap.parse_args(["--help"])

        assert error.type == SystemExit

        captured = capsys.readouterr()
        assert "bar" in captured.out

    @pytest.mark.parametrize(
        'dest, default',
        [("debug", False), ("logfile", None)])
    def test_args(self, ap, dest, default):
        """Test for added run arguments."""
        args = ap.parse_known_args()[0]

        if default is not None:
            assert args.__dict__[dest] == default
        else:
            assert args.__dict__[dest] is None


class Test_convert_analysis_parameters:
    """Test class for converting analysis parameters."""

    @pytest.fixture()
    def universe(self):
        """Univserse fixture."""
        return Universe(TPR, XTC)

    def test_Atomgroup(self, universe):
        """Test AtomGroup conversion."""
        analysis_parameters = {"atomgroup": "all"}
        test_Universe = convert_analysis_parameters(
            analysis_callable=RMSF,
            analysis_parameters=analysis_parameters,
            reference_universe=universe)

        # `None` is returned if no Universe is created
        assert test_Universe is None
        assert analysis_parameters["atomgroup"] == universe.atoms

    def test_zero_atoms(self, universe):
        """Test error if zero atoms are present after conversion."""
        analysis_parameters = {"atomgroup": "name foo"}
        with pytest.raises(ValueError, match="AtomGroup `-atomgroup` with "):
            convert_analysis_parameters(
                analysis_callable=RMSF,
                analysis_parameters=analysis_parameters,
                reference_universe=universe)

    def test_Universe_creation(self, universe):
        """Test the universe creation from different parameters."""
        analysis_parameters = {}
        analysis_parameters["topology_u"] = TPR
        analysis_parameters["topology_format_u"] = None
        analysis_parameters["coordinates_u"] = XTC
        analysis_parameters["trajectory_format_u"] = None
        analysis_parameters["atom_style_u"] = None

        test_u = convert_analysis_parameters(
            analysis_callable=EinsteinMSD,
            analysis_parameters=analysis_parameters,
            reference_universe=None)

        # Extra parameters should be removed
        assert len(analysis_parameters.keys()) == 1
        assert analysis_parameters["u"] is test_u

        # We can not test accros several universes so we just compare
        # the atomic positions
        assert (universe.atoms.positions == test_u.atoms.positions).all()

    def test_only_set_if_key_exists(self, universe):
        """
        Test if keys exist.

        The docparser does not distinguis between mandatory and positional
        arguments. So we only should set the value if it is inside the dict.
        """
        analysis_parameters = {}
        convert_analysis_parameters(
            analysis_callable=RMSF,
            analysis_parameters=analysis_parameters,
            reference_universe=universe)
        assert not analysis_parameters


class Test_run_analsis:
    """Test class for analyze_data."""

    @pytest.fixture()
    def reference_universe_parameters(self):
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
                     reference_universe_parameters,
                     mandatory_parameters,
                     tmpdir):
        """Test with default arguments."""
        with tmpdir.as_cwd():
            run_analsis(
                analysis_callable=InterRDF,
                reference_universe_parameters=reference_universe_parameters,
                mandatory_analysis_parameters=mandatory_parameters)

    def test_optional_paramaters(self,
                                 reference_universe_parameters,
                                 mandatory_parameters,
                                 tmpdir):
        """Test with optional parameters given."""
        opt_params = {"nbins": 100}
        with tmpdir.as_cwd():
            a = run_analsis(
                analysis_callable=InterRDF,
                reference_universe_parameters=reference_universe_parameters,
                mandatory_analysis_parameters=mandatory_parameters,
                optional_analysis_parameters=opt_params)

        assert a.rdf_settings["bins"] == opt_params["nbins"]

    def test_verbose(self,
                     reference_universe_parameters,
                     mandatory_parameters,
                     tmpdir,
                     caplog):
        """Test for being verbose."""
        run_parameters = {"verbose": True}
        caplog.set_level(logging.INFO)
        with tmpdir.as_cwd():
            run_analsis(
                analysis_callable=InterRDF,
                reference_universe_parameters=reference_universe_parameters,
                mandatory_analysis_parameters=mandatory_parameters,
                run_parameters=run_parameters)
        assert "INFO     " in caplog.text

    def test_custom_output(self,
                           reference_universe_parameters,
                           mandatory_parameters,
                           tmpdir):
        """Test for custom output."""
        output_parameters = {"output_directory": "foo",
                             "output_prefix": "bar"}

        with tmpdir.as_cwd():
            os.mkdir("foo")
            run_analsis(
                analysis_callable=InterRDF,
                reference_universe_parameters=reference_universe_parameters,
                mandatory_analysis_parameters=mandatory_parameters,
                output_parameters=output_parameters)

            os.path.isfile(os.path.join("foo",
                                        "bar_InterRDF_count_bins_rdf.csv"))

    def test_output_directory(self,
                              reference_universe_parameters,
                              mandatory_parameters,
                              tmpdir):
        """Test for custom output directory."""
        output_parameters = {"output_directory": "foo"}
        with tmpdir.as_cwd():
            os.mkdir("foo")
            run_analsis(
                analysis_callable=InterRDF,
                reference_universe_parameters=reference_universe_parameters,
                mandatory_analysis_parameters=mandatory_parameters,
                output_parameters=output_parameters)

            os.path.isfile(os.path.join("foo", "InterRDF_count_bins_rdf.csv"))

    def test_output_prefix(self,
                           reference_universe_parameters,
                           mandatory_parameters,
                           tmpdir):
        """Test for custom prefix."""
        output_parameters = {"output_prefix": "foo"}
        with tmpdir.as_cwd():
            os.mkdir("foo")
            run_analsis(
                analysis_callable=InterRDF,
                reference_universe_parameters=reference_universe_parameters,
                mandatory_analysis_parameters=mandatory_parameters,
                output_parameters=output_parameters)

            os.path.isfile("foo_InterRDF_count_bins_rdf.csv")
