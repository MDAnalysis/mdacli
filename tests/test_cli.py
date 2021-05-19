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

from mdacli.cli import setup_clients


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
