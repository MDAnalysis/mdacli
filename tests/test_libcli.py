#!/usr/bin/env python3
#
# Copyright (c) 2020 Authors and contributors
#
# Released under the GNU Public Licence, v2 or any higher version
# SPDX-License-Identifier: GPL-2.0-or-later
"""Test libcli."""

import argparse
import logging
import os
import sys
from io import StringIO
from json.decoder import JSONDecodeError
from pathlib import Path
from unittest.mock import patch

import pytest
from MDAnalysis.analysis import __all__
from MDAnalysis.analysis.base import AnalysisBase
from MDAnalysis.analysis.helix_analysis import HELANAL
from MDAnalysis.analysis.lineardensity import LinearDensity
from MDAnalysis.analysis.msd import EinsteinMSD
from MDAnalysis.analysis.rdf import InterRDF, InterRDF_s
from MDAnalysis.analysis.rms import RMSF
from MDAnalysis.core.universe import Universe
from MDAnalysisTests.core.test_universe import CHOL_GRO
from MDAnalysisTests.datafiles import TPR, XTC
from MDAnalysisTests.topology.test_lammpsdata import LAMMPS_NORESID
from numpy.testing import assert_equal

from mdacli.libcli import (
    KwargsDict,
    add_cli_universe,
    add_output_group,
    add_run_group,
    convert_analysis_parameters,
    create_cli,
    create_universe,
    find_classes_in_modules,
    find_cls_members,
    init_base_argparse,
    run_analysis,
    setup_clients,
    split_argparse_into_groups,
)

from . import example_json
from .test_utils import complete_docstring


@pytest.mark.parametrize(
    ("cmd", "expected"),
    [
        ('-d {"key1":1}', {"key1": 1}),
    ],
)
def test_KwargsDict(cmd, expected):
    """Test dict reading action."""
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "-d",
        action=KwargsDict,
        default=None,
    )
    args = ap.parse_args(cmd.split())
    assert args.d == expected


@pytest.mark.parametrize(
    ("cmd", "expected"),
    [
        (f"-d {os.fspath(example_json)}", {"key1": 1}),
    ],
)
def test_KwargsDict_from_file(cmd, expected):
    """Test dict reading action."""
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "-d",
        action=KwargsDict,
        default=None,
    )
    args = ap.parse_args(cmd.split())
    assert args.d == expected


@pytest.mark.parametrize(
    ("s", "error", "msg"),
    [
        ("-d {fail}", JSONDecodeError, "An error ocurred when reading"),
        ("-d fail", FileNotFoundError, "No such file or directory"),
    ],
)
def test_KwargsDict_error(s, error, msg):
    """Test error."""
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "-d",
        action=KwargsDict,
        default=None,
    )
    with pytest.raises(error, match=msg):
        ap.parse_args(s.split())


def test_find_cls_members():
    """Test several input modules."""
    names = ["helix_analysis", "lineardensity"]
    members = find_cls_members(
        AnalysisBase, [f"MDAnalysis.analysis.{m}" for m in names]
    )
    assert members[0] is HELANAL
    assert members[1] is LinearDensity


@pytest.mark.parametrize(
    "cls", [AnalysisBase, [AnalysisBase, AnalysisBase], (AnalysisBase, AnalysisBase)]
)
def test_find_cls_members_single(cls):
    """Test one input module."""
    members = find_cls_members(cls, ["MDAnalysis.analysis.rdf"])

    assert members[0] is InterRDF
    assert members[1] is InterRDF_s


@pytest.mark.parametrize(
    "cls", [AnalysisBase, [AnalysisBase, AnalysisBase], (AnalysisBase, AnalysisBase)]
)
def test_find_classes_in_modules(cls):
    """Test for finding classes in modules."""
    members = find_classes_in_modules(cls, "MDAnalysis.analysis.rdf")

    assert members[0] is InterRDF
    assert members[1] is InterRDF_s


def test_find_classes_in_modules_None():
    """Test that no module is found."""
    members = find_classes_in_modules(AnalysisBase, "MDAnalysis")

    assert members is None


def test_split_argparse_into_groups():
    """Test splitting argparse Namespace into several dicts."""
    parser = argparse.ArgumentParser()

    group1 = parser.add_argument_group("group1")
    group1.add_argument("--test1", help="test1")

    group2 = parser.add_argument_group("group2")
    group2.add_argument("--test2", help="test2")

    args = parser.parse_args(["--test1", "one", "--test2", "two"])

    arg_grouped_dict = split_argparse_into_groups(parser, args)

    assert arg_grouped_dict["group1"]["test1"] == "one"
    assert arg_grouped_dict["group2"]["test2"] == "two"


@pytest.mark.parametrize(
    ("dest", "default"),
    [("start", None), ("stop", None), ("step", "1"), ("verbose", False)],
)
def test_add_run_group_args(dest, default):
    """Test for added run arguments."""
    parser = argparse.ArgumentParser()
    add_run_group(parser)
    args = parser.parse_args([])
    assert args.__dict__[dest] == default


def test_add_run_group_group(capsys):
    """Test is a run group is added."""
    parser = argparse.ArgumentParser()
    add_run_group(parser)

    parser.print_help()
    captured = capsys.readouterr()

    assert "Analysis Run Parameters:" in captured.out


@pytest.mark.parametrize(
    ("dest", "default"), [("output_prefix", ""), ("output_directory", ".")]
)
def test_add_output_group_arguments(dest, default):
    """Test for added output arguments."""
    parser = argparse.ArgumentParser()
    add_output_group(parser)

    args = parser.parse_known_args()[0]

    if default is not None:
        assert args.__dict__[dest] == default
    else:
        assert args.__dict__[dest] is None


def test_add_output_group_group(capsys):
    """Test is a output group is added."""
    parser = argparse.ArgumentParser()
    add_output_group(parser)

    parser.print_help()
    captured = capsys.readouterr()

    assert "Output Parameters:" in captured.out


@pytest.mark.parametrize("name", ["", "u"])
@pytest.mark.parametrize(
    ("dest", "default"),
    [
        ("topology", "topol.tpr"),
        ("topology_format", None),
        ("coordinates", None),
        ("atom_style", None),
        ("trajectory_format", None),
        ("dimensions", None),
    ],
)
def test_add_cli_universe(name, dest, default):
    """Test for added output arguments."""
    parser = argparse.ArgumentParser()
    add_cli_universe(parser, name)

    args = parser.parse_known_args()[0]

    name = f"_{name}" if name else ""

    if default is not None:
        assert args.__dict__[f"{dest}{name}"] == default
    else:
        assert args.__dict__[f"{dest}{name}"] is None


@pytest.mark.parametrize(
    ("opt", "dest", "val"),
    [
        ("-s", "topology", "foo"),
        ("-top", "topology_format", "foo"),
        ("-f", "coordinates", ["foo", "bar"]),
        ("-traj", "trajectory_format", "foo"),
        ("-atom_style", "atom_style", "foo"),
        ("-dimensions", "dimensions", [10.0, 10.0, 10.0]),
        ("-b", "start", 42),
        ("-e", "stop", 42),
        ("-dt", "step", 42),
    ],
)
def test_setup_clients(opt, dest, val):
    """Test all additional arguments."""
    testargs = ["mdacli", "rmsf", "-atomgroup", "all", opt]
    if type(val) is list:
        for i in val:
            testargs.append(str(i))
    else:
        testargs.append(str(val))

    members = find_cls_members(
        modules=[f"MDAnalysis.analysis.{m}" for m in __all__], cls=AnalysisBase
    )
    actual_mods = ["RMSF", "InterRDF", "LinearDensity"]
    members = [mem for mem in members if mem.__name__ in actual_mods]

    with patch.object(sys, "argv", testargs):
        ap = argparse.ArgumentParser()
        setup_clients(ap, title="title", members=members)
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

    def test_trajectory_format(self):
        """Test non default trajectpry format."""
        u = create_universe(
            topology=StringIO(CHOL_GRO),
            coordinates=StringIO(CHOL_GRO),
            topology_format="GRO",
            trajectory_format="GRO",
        )
        assert len(u.atoms) == 8

    def test_atom_style(self):
        """Test custom atom style."""
        u = create_universe(
            topology=StringIO(LAMMPS_NORESID),
            topology_format="data",
            atom_style="id type x y z",
        )
        assert u.atoms[0].mass == 28.0

    def test_dimensions(self):
        """Test dimensions overwrite."""
        dimensions = [10, 10, 10, 20, 90, 10]
        u = create_universe(TPR, XTC, dimensions=dimensions)
        assert_equal(u.dimensions, dimensions)

    def test_dimensions_no_angles(self):
        """Test dimesions overwrite without giving angles."""
        dimensions = [10, 10, 10]
        u = create_universe(TPR, XTC, dimensions=dimensions)
        assert_equal(u.dimensions, dimensions + [90, 90, 90])

    def test_dimensions_length_error(self):
        """Test dimesions overwrite with wrong number of indices."""
        dimensions = [10, 10]
        with pytest.raises(IndexError, match="at least 3 entries"):
            create_universe(TPR, XTC, dimensions=dimensions)


class Test_init_base_argparse:
    """Test for basic argument parser."""

    @pytest.fixture
    def ap(self):
        """Return the basic parser."""
        return init_base_argparse(name="foo", version="0.0.0", description="bar")

    def test_version(self, ap, capsys):
        """Test version option."""
        with pytest.raises(SystemExit) as error:
            ap.parse_args(["--version"])

        assert error.type is SystemExit

        captured = capsys.readouterr()
        assert captured.out == "foo 0.0.0\n"

    def test_description(self, ap, capsys):
        """Test version option."""
        with pytest.raises(SystemExit) as error:
            ap.parse_args(["--help"])

        assert error.type is SystemExit

        captured = capsys.readouterr()
        assert "bar" in captured.out

    @pytest.mark.parametrize(("dest", "default"), [("debug", False), ("logfile", None)])
    def test_args(self, ap, dest, default):
        """Test for added run arguments."""
        args = ap.parse_known_args()[0]

        if default is not None:
            assert args.__dict__[dest] == default
        else:
            assert args.__dict__[dest] is None


class Test_convert_analysis_parameters:
    """Test class for converting analysis parameters."""

    @pytest.fixture
    def universe(self):
        """Univserse fixture."""
        return Universe(TPR, XTC)

    def test_Atomgroup(self, universe):
        """Test AtomGroup conversion."""
        analysis_parameters = {"atomgroup": "all"}
        test_Universe = convert_analysis_parameters(
            analysis_callable=RMSF,
            analysis_parameters=analysis_parameters,
            reference_universe=universe,
        )

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
                reference_universe=universe,
            )

    def test_None(self, universe):
        """Test that `None` objects are NOT converted.

        Could be default arguments
        """
        analysis_parameters = {"atomgroup": None}
        convert_analysis_parameters(
            analysis_callable=RMSF,
            analysis_parameters=analysis_parameters,
            reference_universe=universe,
        )

        assert analysis_parameters["atomgroup"] is None

    def test_Universe_creation(self, universe):
        """Test the universe creation from different parameters."""
        analysis_parameters = {}
        analysis_parameters["topology_u"] = TPR
        analysis_parameters["topology_format_u"] = None
        analysis_parameters["coordinates_u"] = XTC
        analysis_parameters["trajectory_format_u"] = None
        analysis_parameters["atom_style_u"] = None
        analysis_parameters["dimensions_u"] = None

        test_u = convert_analysis_parameters(
            analysis_callable=EinsteinMSD,
            analysis_parameters=analysis_parameters,
            reference_universe=None,
        )

        # Extra parameters should be removed
        assert len(analysis_parameters.keys()) == 1
        assert analysis_parameters["u"] is test_u

        # We can not test accros several universes so we just compare
        # the atomic positions
        assert (universe.atoms.positions == test_u.atoms.positions).all()

    def test_only_set_if_key_exists(self, universe):
        """Test if keys exist.

        The docparser does not distinguis between mandatory and positional
        arguments. So we only should set the value if it is inside the dict.
        """
        analysis_parameters = {}
        convert_analysis_parameters(
            analysis_callable=RMSF,
            analysis_parameters=analysis_parameters,
            reference_universe=universe,
        )
        assert not analysis_parameters

    def test_multi_atomgroup(self, universe):
        """Test conversion of a list containing multiple AtomGroups."""
        selection_keywords = ["name OW", "name HW*"]
        analysis_parameters = {"p0": ["name OW", "name HW*"]}

        convert_analysis_parameters(
            analysis_callable=complete_docstring,
            analysis_parameters=analysis_parameters,
            reference_universe=universe,
        )

        for i, sel in enumerate(selection_keywords):
            assert universe.select_atoms(sel) == analysis_parameters["p0"][i]

    def test_multi_atomgroup_fail(self, universe):
        """Test error raise for empty atomgroup for multiple AtomGroups."""
        analysis_parameters = {"p0": ["name foo"]}

        with pytest.raises(ValueError, match="AtomGroup `-p0` with "):
            convert_analysis_parameters(
                analysis_callable=complete_docstring,
                analysis_parameters=analysis_parameters,
                reference_universe=universe,
            )

    def test_multi_atomgroup_None(self, universe):
        """Test that `None` objects are NOT converted.

        Could be default arguments
        """
        analysis_parameters = {"p0": None}

        convert_analysis_parameters(
            analysis_callable=complete_docstring,
            analysis_parameters=analysis_parameters,
            reference_universe=universe,
        )

        assert analysis_parameters["p0"] is None


class Test_run_analysis:
    """Test class for analyze_data."""

    @pytest.fixture
    def reference_universe_parameters(self):
        """Universe fixture."""
        kwargs = {}
        kwargs["topology"] = TPR
        kwargs["topology_format"] = None
        kwargs["coordinates"] = XTC
        kwargs["trajectory_format"] = None
        kwargs["atom_style"] = None
        return kwargs

    @pytest.fixture
    def mandatory_parameters(self):
        """Mandatory parameter fixture."""
        kwargs = {}
        kwargs["g1"] = "resid 500"  # water molecules
        kwargs["g2"] = "resid 501"  # water molecules
        return kwargs

    def test_default(self, reference_universe_parameters, mandatory_parameters, tmpdir):
        """Test with default arguments."""
        with tmpdir.as_cwd():
            run_analysis(
                analysis_callable=InterRDF,
                reference_universe_parameters=reference_universe_parameters,
                mandatory_analysis_parameters=mandatory_parameters,
            )

    def test_optional_paramaters(
        self, reference_universe_parameters, mandatory_parameters, tmpdir
    ):
        """Test with optional parameters given."""
        opt_params = {"nbins": 100}
        with tmpdir.as_cwd():
            a = run_analysis(
                analysis_callable=InterRDF,
                reference_universe_parameters=reference_universe_parameters,
                mandatory_analysis_parameters=mandatory_parameters,
                optional_analysis_parameters=opt_params,
            )

        assert a.rdf_settings["bins"] == opt_params["nbins"]

    def test_verbose(
        self, reference_universe_parameters, mandatory_parameters, tmpdir, caplog
    ):
        """Test for being verbose."""
        run_parameters = {"verbose": True}
        caplog.set_level(logging.INFO)
        with tmpdir.as_cwd():
            run_analysis(
                analysis_callable=InterRDF,
                reference_universe_parameters=reference_universe_parameters,
                mandatory_analysis_parameters=mandatory_parameters,
                run_parameters=run_parameters,
            )
        assert "INFO     " in caplog.text

    def test_custom_output(
        self, reference_universe_parameters, mandatory_parameters, tmpdir
    ):
        """Test for custom output."""
        output_parameters = {"output_directory": "foo", "output_prefix": "bar"}

        with tmpdir.as_cwd():
            Path.mkdir("foo")
            run_analysis(
                analysis_callable=InterRDF,
                reference_universe_parameters=reference_universe_parameters,
                mandatory_analysis_parameters=mandatory_parameters,
                output_parameters=output_parameters,
            )

            assert (Path("foo") / "bar_InterRDF_count_bins_rdf.csv").is_file()

    def test_output_directory(
        self, reference_universe_parameters, mandatory_parameters, tmpdir
    ):
        """Test for custom output directory."""
        output_parameters = {"output_directory": "foo"}
        with tmpdir.as_cwd():
            Path.mkdir("foo")
            run_analysis(
                analysis_callable=InterRDF,
                reference_universe_parameters=reference_universe_parameters,
                mandatory_analysis_parameters=mandatory_parameters,
                output_parameters=output_parameters,
            )

            assert (Path("foo") / "InterRDF_count_bins_rdf.csv").is_file()

    def test_output_prefix(
        self, reference_universe_parameters, mandatory_parameters, tmpdir
    ):
        """Test for custom prefix."""
        output_parameters = {"output_prefix": "foo"}
        with tmpdir.as_cwd():
            Path.mkdir("foo")
            run_analysis(
                analysis_callable=InterRDF,
                reference_universe_parameters=reference_universe_parameters,
                mandatory_analysis_parameters=mandatory_parameters,
                output_parameters=output_parameters,
            )

            assert Path("foo_InterRDF_count_bins_rdf.csv").is_file()


class Test_create_cli:
    """Test class for CLI creation."""

    @pytest.fixture
    def parameters(self):
        """Inject key,value pair into parameter dictionary."""
        return {
            "callable": lambda x: x,
            "desc": "foo",
            "desc_long": "bar",
            "positional": {"pp": {"desc": "pp desc"}},
            "optional": {"po": {"desc": "po desc"}},
        }

    def cli(self, parameters):
        """Inject argument into subparser."""
        ap = argparse.ArgumentParser()
        subparser = ap.add_subparsers()

        create_cli(sub_parser=subparser, interface_name="foo", parameters=parameters)
        return subparser.choices["foo"]

    def test_description(self, parameters):
        """Test parser description."""
        cli = self.cli(parameters)

        desc = parameters["desc"]
        desc += "\n\n" + parameters["desc_long"]
        assert cli.description == desc

    @pytest.mark.parametrize(
        "attr",
        ["start", "stop", "step", "verbose", "output_prefix", "output_directory"],
    )
    def test_common_args(self, parameters, attr):
        """Test common cli parameters."""
        cli = self.cli(parameters)

        args = cli.parse_known_args(["-pp", "foo"])[0]

        getattr(args, attr)

    def test_add_output_group(self, parameters):
        """Test if `save` not added."""
        callable = type("", (), {})()
        callable.save = True
        parameters["callable"] = callable

        cli = self.cli(parameters)
        args = cli.parse_known_args(["-pp", "foo"])[0]

        with pytest.raises(AttributeError):
            _ = args.output_directory

    @pytest.mark.parametrize("argument", ["positional", "optional"])
    @pytest.mark.parametrize(
        ("val_type", "arg_type", "value"),
        [
            ("bool", None, False),
            ("str", str, "foo"),
            ("list", list, None),
            ("tuple", tuple, None),
            ("dict", None, None),
            ("int", int, 42),
            ("float", float, 1.1),
            ("complex", complex, 1j),
            ("NoneType", None, None),
            ("AtomGroup", str, None),
            ("MDAnalysis.core.groups.AtomGroup", str, None),
            ("list[AtomGroup]", str, None),
        ],
    )
    def test_arguments(self, parameters, argument, val_type, arg_type, value):
        """Test for existance and default value of arguments."""
        help = "p0 desc"
        opt_params = {"p0": {"type": val_type, "desc": help}}
        if argument == "optional":
            opt_params["p0"]["default"] = value
        parameters[argument] = opt_params

        cli = self.cli(parameters)
        # last or pre last element in list is p0
        action = cli._actions[-1] if argument == "optional" else cli._actions[-2]

        assert action.dest == "p0"
        assert action.option_strings[0] == "-p0"
        assert help in action.help
        assert action.type == arg_type
        if argument == "optional":
            assert action.default == value
        elif argument == "positional":
            assert action.required

    def test_true_bool_argument(self, parameters):
        """Test if no is added for bool that if swapped from true to false."""
        opt_params = {"p0": {"type": "bool", "desc": "p0 desc", "default": True}}
        parameters["optional"] = opt_params

        cli = self.cli(parameters)
        action = cli._actions[-1]

        assert action.option_strings[0] == "-no-p0"

    @pytest.mark.xfail(reason="TODO check choices parsing for numbers")
    @pytest.mark.parametrize(
        ("val", "type", "choices"),
        [
            ("{'1','2','-1'},", int, [1, 2, -1]),
            ("{'1E8','-1.4'},", float, [1e8, -1.4]),
            ("{ 'a', 'b','c' },", str, ["a", "b", "c"]),
        ],
    )
    def tests_choices(self, parameters, val, type, choices):
        """Test choices get parsed correctly."""
        opt_params = {"p0": {"type": val, "desc": "p0 desc", "default": True}}
        parameters["optional"] = opt_params

        cli = self.cli(parameters)
        action = cli._actions[-1]
        assert action.choices == choices
        assert action.type is type

    def test_atomgroup_extra_argument(self, parameters):
        """Test if a universe group is added."""
        opt_params = {"p0": {"type": "AtomGroup", "desc": "p0 desc", "default": None}}
        parameters["optional"] = opt_params

        cli = self.cli(parameters)

        titles = [a.title for a in cli._action_groups]
        assert "Reference Universe Parameters" in titles

    def test_universe_argument(self, parameters):
        """Test if universe arguments are added."""
        opt_params = {"p0": {"type": "Universe", "desc": "p0 desc", "default": None}}
        parameters["optional"] = opt_params

        cli = self.cli(parameters)
        args = cli.parse_known_args(["-pp", "foo"])[0]

        _ = args.topology_p0

    @pytest.mark.parametrize("val_type", ["list", "tuple", "list[AtomGroup]"])
    def test_iterable(self, parameters, val_type):
        """Test if iterable types has nargs='+' attribute."""
        opt_params = {"p0": {"type": val_type, "desc": ""}}
        parameters["positional"] = opt_params

        cli = self.cli(parameters)
        action = cli._actions[-2]

        assert action.nargs == "+"
