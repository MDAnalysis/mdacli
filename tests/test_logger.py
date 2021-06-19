#!/usr/bin/env python3
# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
#
# Copyright (c) 2021 Authors and contributors
#
# Released under the GNU Public Licence, v2 or any higher version
# SPDX-License-Identifier: GPL-2.0-or-later
"""Test mdacli logger."""
import logging

import mdacli.logger


class Test_setup_logger:
    """Test the setup_logger."""

    def log_msg(self):
        """Info message parsed to the log."""
        logger = logging.getLogger("test")
        logger.info("foo")

    def test_default_log(self, caplog):
        """Default message in STDOUT."""
        caplog.set_level(logging.INFO)
        with mdacli.logger.setup_logging(logfile=None, debug=False):
            self.log_msg()
            assert "foo" in caplog.text

    def test_info_log(self, tmpdir, caplog):
        """Default message in STDOUT and file."""
        caplog.set_level(logging.INFO)
        with tmpdir.as_cwd():
            with mdacli.logger.setup_logging(logfile="logfile", debug=False):
                self.log_msg()
                assert "foo" in caplog.text
            with open("logfile.log", "r") as f:
                log = f.read()

            assert log in caplog.text

    def test_debug_log(self, tmpdir, caplog):
        """Debug message in STDOUT and file."""
        caplog.set_level(logging.INFO)
        with tmpdir.as_cwd():
            with mdacli.logger.setup_logging(logfile="logfile", debug=True):
                self.log_msg()
                assert "test:test_logger.py:20 foo\n" in caplog.text
                with open("logfile.log", "r") as f:
                    log = f.read()

                assert "test_logger.py:test:log_msg:20:foo\n" in log
