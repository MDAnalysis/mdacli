#!/usr/bin/env python
#
# Copyright (c) 2024 Authors and contributors
# (see the AUTHORS.rst file for the full list of names)
#
# Released under the GNU Public Licence, v3 or any higher version
# SPDX-License-Identifier: GPL-3.0-or-later
"""Test module for mdacli."""

from MDAnalysis.analysis.base import AnalysisBase

from mdacli import cli


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
