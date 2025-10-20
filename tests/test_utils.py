#!/usr/bin/env python3
#
# Copyright (c) 2021 Authors and contributors
#
# Released under the GNU Public Licence, v2 or any higher version
# SPDX-License-Identifier: GPL-2.0-or-later
"""Test mdacli utils."""

from collections import defaultdict
from math import isclose

import pytest

from mdacli.utils import (
    _exit_if_a_is_b,
    convert_str_time,
    parse_callable_signature,
    parse_docs,
    split_time_unit,
)


def test__exit_if_a_is_b():
    """Test for a SystemExit using pytest.raises."""
    msg = "foo"
    with pytest.raises(SystemExit, match=msg) as error:
        _exit_if_a_is_b(1, 1, msg=msg)
    assert error.type is SystemExit


@pytest.mark.parametrize(
    ("x", "frame"),
    [
        ("1", 1),
        ("-1", -1),
        ("1ns", 1000),
        ("1.5ns", 1500),
        ("12e3", 12e3),
        ("12e3ps", 12e3),
    ],
)
def test_convert_str_time(x, frame):
    """Test convert string to time."""
    assert frame == convert_str_time(x, dt=1)


def test_convert_str_time_dt():
    """Test convert string to time in ps."""
    assert convert_str_time("10ps", dt=10) == 1


def test_convert_str_time_raise():
    """Test convert string to time ValueError."""
    with pytest.raises(ValueError, match="Only integers or time step combinations"):
        convert_str_time("0.1", dt=1)


@pytest.mark.parametrize(
    ("s", "expected"),
    [
        ("00", (0.0, "")),
        ("0", (0.0, "")),
        ("-0", (0.0, "")),
        ("0.", (0.0, "")),
        (".0ns", (0.0, "ns")),
        ("12.3ps", (12.3, "ps")),
        ("12.34minute", (12.34, "minute")),
        ("-12.34", (-12.34, "")),
        ("10E40sd", (float("10E40"), "sd")),
        ("1E4", (float("1E4"), "")),
        ("1E-4hours", (float("1E-4"), "hours")),
        ("-10E40nanosecond", (float("-10E40"), "nanosecond")),
        ("-10E-40", (float("-10E-40"), "")),
        ("-.1E-4", (float("-.1E-4"), "")),
        ("10.2E30", (float("10.2E30"), "")),
        (".10E30", (float(".10E30"), "")),
        ("-10.2E30", (float("-10.2E30"), "")),
        ("-.10E30", (float("-.10E30"), "")),
    ],
)
def test_string_to_timestep(s, expected):
    """Test correct conversion."""
    time_, unit_ = split_time_unit(s)
    assert isclose(time_, expected[0])
    assert unit_ == expected[1]


@pytest.mark.parametrize(
    "s",
    [
        "12.34 # with comment",
        "12.34.12fdsfdsfds",
        "1E4.4",
        ".10E30E",
        "10.2.E30",
        "123#with wrong comment",
        "E10",
        ".E10",
        "-e",
        "-E10",
        "10-10E19",
    ],
)
def test_string_to_timestep_wrong(s):
    """Test correct conversion."""
    with pytest.raises(IndexError):
        split_time_unit(s)


def complete_docstring(p0, p1="foo", p2=True, p3=42):  # NOQA: ARG001
    """One-line description.

    Multi-
    line-
    description.

    Parameters
    ----------
    p0 : list[AtomGroup]
        Param 0 is a list
    p1 : str or int
        Param 1 description.
    p2 : bool
        Param 2
        description.
    p3 : int
        Param 3
        description.
    """
    summary = "One-line description."
    summary_extended = "Multi-\nline-\ndescription."
    params = defaultdict(
        dict,
        {
            "p3": {"type": "int", "desc": "Param 3 description."},
            "p2": {"type": "bool", "desc": "Param 2 description."},
            "p1": {"type": "str", "desc": "Param 1 description."},
            "p0": {"type": "list[AtomGroup]", "desc": "Param 0 is a list"},
        },
    )
    return summary, summary_extended, params


def no_long_docstring(p0, p1="foo", p2=True):  # NOQA: ARG001
    """One-line description.

    Parameters
    ----------
    p0 : list[AtomGroup]
        Param 0 is a list
    p1 : str or int
        Param 1 description.
    p2 : bool
        Param 2
        description.
    """
    summary = "One-line description."
    summary_extended = ""
    params = defaultdict(
        dict,
        {
            "p2": {"type": "bool", "desc": "Param 2 description."},
            "p1": {"type": "str", "desc": "Param 1 description."},
            "p0": {"type": "list[AtomGroup]", "desc": "Param 0 is a list"},
        },
    )
    return summary, summary_extended, params


@pytest.mark.parametrize("func", [complete_docstring, no_long_docstring])
def test_parse_docstring(func):
    """Test doc string parsing."""
    summary, summary_extended, params = func(1)
    summary_parse, summary_extended_parse, params_parse = parse_docs(func)

    assert summary == summary_parse
    assert summary_extended == summary_extended_parse
    assert params == params_parse


def test_parse_callable_signature():
    """Test callable signature parsing."""
    parameters = parse_callable_signature(complete_docstring)
    summary, summary_extended, params = parse_docs(complete_docstring)
    optional = {"p1": params["p1"], "p2": params["p2"], "p3": params["p3"]}
    optional["p1"]["default"] = "foo"
    optional["p2"]["default"] = True
    optional["p3"]["default"] = 42

    assert parameters["callable"] == complete_docstring
    assert parameters["positional"] == {"p0": params["p0"]}
    assert parameters["optional"] == optional
    assert parameters["desc"] == summary
    assert parameters["desc_long"] == summary_extended
