#!/usr/bin/env python3
# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
#
# Copyright (c) 2024 Authors and contributors
#
# Released under the GNU Public Licence, v2 or any higher version
# SPDX-License-Identifier: GPL-2.0-or-later
"""Mock module for mdacli to test logging."""

import logging

from MDAnalysis.analysis.base import AnalysisBase


logger = logging.getLogger(__name__)


class Tester(AnalysisBase):
    """Mock class for mdacli. Implements only the minimum requirements.

    Currently only logs messages at different levels to check the verbosity and
    debug flags in the CLI.

    Parameters
    ----------
    atomgroup : AtomGroup or Universe
    """

    def __init__(self, atomgroup, **kwargs):
        """Initialise the Tester class."""
        super(Tester, self).__init__(atomgroup.universe.trajectory, **kwargs)
        logger.info("This is an info message")
        logger.warn("This is a warning")
        logger.debug("This is a debug message")

    def _prepare(self):
        """Prepare the analysis."""
        pass

    def _single_frame(self):
        """Analyse a single frame."""
        pass

    def _conclude(self):
        """Conclude the analysis."""
        pass
