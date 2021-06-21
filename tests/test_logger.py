#!/usr/bin/env python3
# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
#
# Copyright (c) 2021 Authors and contributors
#
# Released under the GNU Public Licence, v2 or any higher version
# SPDX-License-Identifier: GPL-2.0-or-later
"""Test mdacli logger."""
import logging

import pytest

import mdacli.logger


class Test_setup_logger:
    """Test the setup_logger."""

    def test_default_log(self, caplog):
        """Default message in STDOUT."""
        caplog.set_level(logging.INFO)
        logger = logging.getLogger("test")
        with mdacli.logger.setup_logging(logger,
                                         logfile=None,
                                         debug=False):
            logger.info("foo")
            assert "foo" in caplog.text

    def test_info_log(self, tmpdir, caplog):
        """Default message in STDOUT and file."""
        caplog.set_level(logging.INFO)
        logger = logging.getLogger("test")
        with tmpdir.as_cwd():
            # Explicityly leave out the .dat file ending to check if this
            # is created by the function.
            with mdacli.logger.setup_logging(logger,
                                             logfile="logfile",
                                             debug=False):
                logger.info("foo")
                assert "foo" in caplog.text
            with open("logfile.log", "r") as f:
                log = f.read()

            assert "foo" in log

    def test_debug_log(self, tmpdir, caplog):
        """Debug message in STDOUT and file."""
        caplog.set_level(logging.INFO)
        logger = logging.getLogger("test")
        with tmpdir.as_cwd():
            with mdacli.logger.setup_logging(logger,
                                             logfile="logfile",
                                             debug=True):
                logger.info("foo")
                assert "test:test_logger.py:54 foo\n" in caplog.text
                
            with open("logfile.log", "r") as f:
                log = f.read()

            assert "test_logger.py:test:test_debug_log:54: foo\n" in log
