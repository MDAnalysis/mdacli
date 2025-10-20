#!/usr/bin/env python3
#
# Copyright (c) 2021 Authors and contributors
#
# Released under the GNU Public Licence, v2 or any higher version
# SPDX-License-Identifier: GPL-2.0-or-later
"""Command line interfaces for MDAnalysis analysis classes."""

from mdacli.cli import cli

from ._version import __version__  # noqa: F401

__all__ = ["cli"]
__authors__ = "MDAnalysis Development Team and contributors"
