#!/usr/bin/env python3
# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
#
# Copyright (c) 2021 Authors and contributors
#
# Released under the GNU Public Licence, v2 or any higher version
# SPDX-License-Identifier: GPL-2.0-or-later
"""Test mdacli utils."""
import argparse
import os
import pytest
from json import JSONDecodeError

from mdacli.cli import convert_str_time
from mdacli import libcli

from . import example_json

@pytest.mark.parametrize('x, frame',
                         (('1', 1),
                          ('-1', -1),
                          ('1ns', 1000),
                          ('12e3', 12e3),
                          ('12e3ps', 12e3)))
def test_convert_str_time(x, frame):
    """Test convert string to time."""
    assert frame == convert_str_time(x, dt=1)


def test_convert_str_time_dt():
    """Test convert string to time in ps."""
    assert 1 == convert_str_time("10ps", dt=10)


def test_convert_str_time_raise():
    """Test convert string to time ValueError."""
    with pytest.raises(ValueError):
        convert_str_time('0.1', dt=1)


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


@pytest.mark.parametrize(
    's,error',
    [
        ("-d {fail}", JSONDecodeError),
        ("-d fail", FileNotFoundError),
        ]
    )
def test_KwargsDict_error(s, error):
    """Test error."""
    ap = argparse.ArgumentParser()
    ap.add_argument(
        '-d',
        action=libcli.KwargsDict,
        default=None,
        )
    with pytest.raises(error):
        ap.parse_args(s.split())
