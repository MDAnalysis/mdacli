#!/usr/bin/env python3
# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
#
# Copyright (c) 2021 Authors and contributors
#
# Released under the GNU Public Licence, v2 or any higher version
# SPDX-License-Identifier: GPL-2.0-or-later
"""Test mdacli utils."""
from math import isclose

import pytest

from mdacli.utils import _exit_if_a_is_b, convert_str_time, split_time_unit


def test__exit_if_a_is_b():
    """Test for a SystemExit using pytest.raises."""
    msg = "foo"
    with pytest.raises(SystemExit, match=msg) as error:
        _exit_if_a_is_b(1, 1, msg=msg)
    assert error.type == SystemExit


@pytest.mark.parametrize('x, frame',
                         (('1', 1),
                          ('-1', -1),
                          ('1ns', 1000),
                          ('1.5ns', 1500),
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
    's,expected',
    [
        ['00', (0.0, '')],
        ['0', (0.0, '')],
        ['-0', (0.0, '')],
        ['0.', (0.0, '')],
        ['.0ns', (0.0, 'ns')],
        ['12.3ps', (12.3, 'ps')],
        ['12.34minute', (12.34, 'minute')],
        ['-12.34', (-12.34, '')],
        ['10E40sd', (float('10E40'), 'sd')],
        ['1E4', (float('1E4'), '')],
        ['1E-4hours', (float('1E-4'), 'hours')],
        ['-10E40nanosecond', (float('-10E40'), 'nanosecond')],
        ['-10E-40', (float('-10E-40'), '')],
        ['-.1E-4', (float('-.1E-4'), '')],
        ['10.2E30', (float('10.2E30'), '')],
        ['.10E30', (float('.10E30'), '')],
        ['-10.2E30', (float('-10.2E30'), '')],
        ['-.10E30', (float('-.10E30'), '')],
        ]
    )
def test_string_to_timestep(s, expected):
    """Test correct conversion."""
    time_, unit_ = split_time_unit(s)
    assert isclose(time_, expected[0])
    assert unit_ == expected[1]


@pytest.mark.parametrize(
    's',
    [
        '12.34 # with comment',
        '12.34.12fdsfdsfds',
        '1E4.4',
        '.10E30E',
        '10.2.E30',
        '123#with wrong comment',
        'E10',
        '.E10',
        '-e',
        '-E10',
        '10-10E19',
        ]
    )
def test_string_to_timestep_wrong(s):
    """Test correct conversion."""
    with pytest.raises(IndexError):
        split_time_unit(s)
