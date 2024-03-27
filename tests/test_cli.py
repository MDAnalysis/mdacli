#!/usr/bin/env python3
# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
#
# Copyright (c) 2021 Authors and contributors
#
# Released under the GNU Public Licence, v2 or any higher version
# SPDX-License-Identifier: GPL-2.0-or-later
"""Test mdacli cli."""

import subprocess

import pytest
from MDAnalysisTests.datafiles import TPR, XTC


def test_required_args():
    """Test that there is a module given."""
    with pytest.raises(subprocess.CalledProcessError):
        subprocess.check_call(["mda"])


def test_wrong_module():
    """Test for a non existent module."""
    with pytest.raises(subprocess.CalledProcessError):
        subprocess.check_call(["mda", "foo"])


@pytest.mark.parametrize("args", ("version", "debug", "help"))
def test_extra_options(args):
    """Test for a ab extra option."""
    subprocess.check_call(["mda", "--" + args])


@pytest.mark.parametrize("args", ("RMSF", "rmsf"))
def test_case_insensitive(args):
    """Test for being case insensitive."""
    subprocess.check_call(["mda", args, "-h"])


@pytest.mark.parametrize('args', ("RMSF", "rmsf"))
def test_case_insensitive_with_flags(args):
    """Test for module name being case insensitive with additional flags."""
    # Check if it still works if the module name is not the second argument
    subprocess.check_call(['mda', '--debug', args, "-h"])


def test_running_analysis(tmpdir):
    """Test running a complete analysis."""
    with tmpdir.as_cwd():
        subprocess.check_call(
            ["mda", "rmsf", "-s", TPR, "-f", XTC, "-atomgroup", "all"]
        )


def test_verbosity_level_warning(caplog):
    """Test the log level warning."""
    # This should only print warning messages
    output = subprocess.check_output(
        ["./tests/run_tester",
         "tester", "-s", TPR, "-f", XTC, "-atomgroup", "all"],
        text=True,
    )
    assert "This is a warning" in output
    # Cross-check that info and debug messages are not printed
    assert "This is a debug message" not in output
    assert "This is an info message" not in output


def test_verbosity_level_info(caplog):
    """Test the log level info."""
    # This should only print warning and info messages
    output = subprocess.check_output(
        [
            "./tests/run_tester",
            "tester", "-s", TPR, "-f", XTC,
            "-atomgroup", "all",
            "-v",
        ],
        text=True,
    )
    assert "This is an info message" in output
    assert "This is a warning" in output
    # Cross-check that debug messages are not printed
    assert "This is a debug message" not in output


def test_verbosity_level_debug(caplog):
    """Test the log level debug."""
    # This should print all messages
    output = subprocess.check_output(
        [
            "./tests/run_tester", "--debug",
            "tester", "-s", TPR, "-f", XTC,
            "-atomgroup", "all",
            "-v",
        ],
        text=True,
    )
    assert "This is an info message" in output
    assert "This is a warning" in output
    assert "This is a debug message" in output
