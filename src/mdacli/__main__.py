#!/usr/bin/env python3
# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
#
# Copyright (c) 2021 Authors and contributors
#
# Released under the GNU Public Licence, v2 or any higher version
# SPDX-License-Identifier: GPL-2.0-or-later

"""
A command line interface (CLI) to the analysis modules of MDAnalysis.

The modules are all structured as part of a single mdacli wrapper, and invoked
with commands like `mda RMSD`. This command uses the class
:class:`MDAnalysis.analysis.rms.RMSD` for calculating the RMSD.
Documentation for each module can be found at the respective sections on the
`MDAnalysis Documentation`_, as well as
`mda command -h`.

.. _`MDAnalysis Documentation`:
   https://docs.mdanalysis.org/stable/documentation_pages/analysis_modules.html
"""

from MDAnalysis.analysis import __all__

import mdacli


def main():
    """Execute main CLI entry point."""
    # Thse module are currently not supported. Either due a different
    # structure or due parameters that aree not supported by our parser.
    skip_mods = ['AnalysisFromFunction',
                 'HydrogenBondAnalysis',
                 'WaterBridgeAnalysis',
                 'Contacts',
                 'Dihedral',
                 'PersistenceLength',
                 'InterRDF_s']

    mdacli.cli(name="MDAnalysis",
               module_list=__all__,
               version=mdacli.__version__,
               description=__doc__,
               skip_modules=skip_mods,
               ignore_warnings=True)


if __name__ == '__main__':
    main()
