#!/usr/bin/env python3
#
# Copyright (c) 2021 Authors and contributors
#
# Released under the GNU Public Licence, v2 or any higher version
# SPDX-License-Identifier: GPL-2.0-or-later
"""Test mdacli save."""

import sys
from unittest.mock import patch

import mdacli.save


def test_get_cli_input():
    """Test the cli input format."""
    testargs = ["maicos", "foo", "foo bar"]
    with patch.object(sys, "argv", testargs):
        cli_in = 'Command line was: maicos foo "foo bar"'
        assert mdacli.save.get_cli_input() == cli_in
