#!/usr/bin/env python3
# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
#
# Copyright (c) 2020 Authors and contributors
#
# Released under the GNU Public Licence, v2 or any higher version
# SPDX-License-Identifier: GPL-2.0-or-later
"""Test mdacli."""

import pytest
from MDAnalysis.analysis.rdf import InterRDF
from MDAnalysis.tests.datafiles import DCD, PSF

from mdacli.cli import analyze_data, convert_str_time


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


class Test_analyze_data(object):
    """Test class for analyze_data."""

    @pytest.fixture()
    def kwargs(self):
        """Keyword arguments for run."""
        kwargs = {}
        kwargs["begin"] = "0"
        kwargs["end"] = "1"
        kwargs["dt"] = "1"
        kwargs["verbose"] = False
        kwargs["g1"] = "all"
        kwargs["g2"] = "all"
        kwargs["func"] = None
        return kwargs

    def test_run(self, kwargs, tmpdir):
        """Simple test run."""
        with tmpdir.as_cwd():
            analyze_data(PSF, DCD, analysis_callable=InterRDF, **kwargs)
