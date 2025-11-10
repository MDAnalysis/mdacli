#!/usr/bin/env python3
#
# Copyright (c) 2021 Authors and contributors
#
# Released under the GNU Public Licence, v2 or any higher version
# SPDX-License-Identifier: GPL-2.0-or-later
"""Test mdacli cli."""

import subprocess
import sys
from pathlib import Path

import pytest
from MDAnalysisTests.datafiles import TPR, XTC

tester_class = (Path().absolute() / "tests/run_tester").as_posix()


def test_required_args():
    """Test that there is a module given."""
    with pytest.raises(subprocess.CalledProcessError):
        subprocess.check_call(["mda"])


def test_wrong_module():
    """Test for a non existent module."""
    with pytest.raises(subprocess.CalledProcessError):
        subprocess.check_call(["mda", "foo"])


@pytest.mark.parametrize("args", ["version", "debug", "help"])
def test_extra_options(args):
    """Test for a ab extra option."""
    subprocess.check_call(["mda", "--" + args])


@pytest.mark.parametrize("args", ["RMSF", "rmsf"])
def test_case_insensitive(args):
    """Test for being case insensitive."""
    subprocess.check_call(["mda", args, "-h"])


@pytest.mark.parametrize("args", ["RMSF", "rmsf"])
def test_case_insensitive_with_flags(args):
    """Test for module name being case insensitive with additional flags."""
    # Check if it still works if the module name is not the second argument
    subprocess.check_call(["mda", "--debug", args, "-h"])


def test_subparser_setup_for_tab_completion():
    """Test that subparsers are correctly set up for tab-completion.

    This verifies that RMSF and RMSD modules are registered as subcommands,
    which is what argcomplete needs for tab-completion to work.
    """
    from MDAnalysis.analysis.base import AnalysisBase

    from src.mdacli.libcli import find_cls_members, init_base_argparse, setup_clients

    modules = find_cls_members(AnalysisBase, ["MDAnalysis.analysis.rms"])

    parser = init_base_argparse(
        name="MDAnalysis", version="0.1.0", description="Test CLI"
    )

    setup_clients(parser, title="MDAnalysis Analysis Modules", members=modules)

    subparser_action = [
        a for a in parser._subparsers._group_actions if hasattr(a, "choices")
    ][0]

    choices = list(subparser_action.choices.keys())
    assert "RMSF" in choices
    assert "RMSD" in choices


def test_argcomplete_working():
    """Test that argcomplete is properly registered and working"""
    import argparse
    from unittest.mock import patch

    from MDAnalysis.analysis import __all__

    import src.mdacli
    from src.mdacli.cli import cli

    skip_mods = [
        "AnalysisFromFunction",
        "HydrogenBondAnalysis",
        "WaterBridgeAnalysis",
        "Contacts",
        "PersistenceLength",
        "InterRDF_s",
    ]

    with patch("src.mdacli.cli.argcomplete.autocomplete") as mock_autocomplete:
        # Mock sys.argv to prevent the CLI from actually running
        with patch("sys.argv", ["mda", "--help"]):
            try:
                cli(
                    name="MDAnalysis",
                    module_list=[f"MDAnalysis.analysis.{m}" for m in __all__],
                    version=src.mdacli.__version__,
                    description="Test",
                    skip_modules=skip_mods,
                    ignore_warnings=True,
                )
            except SystemExit:
                pass

    # Verify that argcomplete.autocomplete was called
    assert mock_autocomplete.called, "argcomplete.autocomplete() was not called"

    # Verify it was called with an ArgumentParser instance
    call_args = mock_autocomplete.call_args
    assert call_args is not None, (
        "argcomplete.autocomplete() was called with no arguments"
    )

    parser_arg = call_args[0][0]  # First positional argument
    assert isinstance(parser_arg, argparse.ArgumentParser), (
        f"argcomplete.autocomplete() should be called with ArgumentParser, got {type(parser_arg)}"
    )


def test_running_analysis(tmpdir):
    """Test running a complete analysis."""
    with tmpdir.as_cwd():
        subprocess.check_call(
            ["mda", "rmsf", "-s", TPR, "-f", XTC, "-atomgroup", "all"]
        )


def test_verbosity_level_warning():
    """Test the log level warning."""
    # This should only print warning messages
    output = subprocess.check_output(
        [
            sys.executable,
            tester_class,
            "tester",
            "-s",
            TPR,
            "-f",
            XTC,
            "-atomgroup",
            "all",
        ],
        text=True,
    )
    assert "This is a warning" in output
    # Cross-check that info and debug messages are not printed
    assert "This is a debug message" not in output
    assert "This is an info message" not in output


def test_verbosity_level_info():
    """Test the log level info."""
    # This should only print warning and info messages
    output = subprocess.check_output(
        [
            sys.executable,
            tester_class,
            "tester",
            "-s",
            TPR,
            "-f",
            XTC,
            "-atomgroup",
            "all",
            "-v",
        ],
        text=True,
    )
    assert "This is an info message" in output
    assert "This is a warning" in output
    # Cross-check that debug messages are not printed
    assert "This is a debug message" not in output


def test_verbosity_level_debug():
    """Test the log level debug."""
    # This should print all messages
    output = subprocess.check_output(
        [
            sys.executable,
            tester_class,
            "--debug",
            "tester",
            "-s",
            TPR,
            "-f",
            XTC,
            "-atomgroup",
            "all",
            "-v",
        ],
        text=True,
    )
    assert "This is an info message" in output
    assert "This is a warning" in output
    assert "This is a debug message" in output
