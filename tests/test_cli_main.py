#!/usr/bin/env python3
# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
#
# Copyright (c) 2020 Authors and contributors
#
# Released under the GNU Public Licence, v2 or any higher version
# SPDX-License-Identifier: GPL-2.0-or-later
"""Test mdacli."""
import sys

import pytest


from MDAnalysis.analysis.rdf import InterRDF
from MDAnalysis.tests.datafiles import PSF, DCD

# Workaround since we have no real module
sys.path.append("..")

from mdacli.cli import convert_str_time  # noqa: E402


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


class TestMain(object):

    @pytest.fixture()
    def kwargs(self):
        kwargs = {}
        kwargs["begin"] = "0"
        kwargs["end"] = "1"
        kwargs["dt"] = "1"
        kwargs["verbose"] = False
        kwargs["g1"] = "all"
        kwargs["g2"] = "all"
        kwargs["func"] = None
        return kwargs

    def test_pickling(self, kwargs, tmpdir):
        with tmpdir.as_cwd():
            with pytest.warns(UserWarning):
                main(PSF, DCD, analysis_callable=InterRDF, **kwargs)
                open("{}.pickle".format(InterRDF.__name__))
