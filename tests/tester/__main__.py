#!/usr/bin/env python
# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
#
# Copyright (c) 2024 Authors and contributors
# (see the AUTHORS.rst file for the full list of names)
#
# Released under the GNU Public Licence, v3 or any higher version
# SPDX-License-Identifier: GPL-3.0-or-later
"""Analyse molecular dynamics simulation of interfacial and confined systems."""

from mdacli import cli

from MDAnalysis.analysis.base import AnalysisBase

def main():
    """Execute main CLI entry point."""
    cli(
        name="Tester",
        module_list=["tester"],
        base_class=AnalysisBase,
        version=0.0,
        description="test",
        ignore_warnings=True,
    )


if __name__ == "__main__":
    main()
