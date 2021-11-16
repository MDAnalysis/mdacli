#!/usr/bin/env python3
# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
#
# Copyright (c) 2021 Authors and contributors
#
# Released under the GNU Public Licence, v2 or any higher version
# SPDX-License-Identifier: GPL-2.0-or-later
"""Test mdacli utils."""
import pytest

from mdacli.utils import _exit_if_a_is_b, convert_str_time


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
