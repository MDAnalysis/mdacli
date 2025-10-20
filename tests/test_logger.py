#!/usr/bin/env python3
#
# Copyright (c) 2024 Authors and contributors
#
# Released under the GNU Public Licence, v2 or any higher version
# SPDX-License-Identifier: GPL-2.0-or-later
"""Test mdacli logger."""

import logging
import warnings
from pathlib import Path

import pytest

from mdacli.logger import check_suffix, setup_logging


@pytest.mark.parametrize("filename", ["example.txt", Path("example.txt")])
def test_check_suffix(filename):
    """Check suffix tetsing."""
    result = check_suffix(filename, ".txt")

    assert str(result) == "example.txt"
    assert isinstance(result, type(filename))


@pytest.mark.parametrize("filename", ["example", Path("example")])
def test_warning_on_missing_suffix(filename):
    """Check issued warning in missing suffix."""
    match = r"The file name should have a '\.txt' extension."
    with pytest.warns(UserWarning, match=match):
        result = check_suffix(filename, ".txt")

    assert str(result) == "example.txt"
    assert isinstance(result, type(filename))


def test_warnings_in_log(caplog):
    """Test that warnings are forwarded to the logger.

    Keep this test at the top since it seems otherwise there are some pytest
    issues...
    """
    logger = logging.getLogger()

    with setup_logging(logger):
        warnings.warn("A warning", stacklevel=1)

    assert "A warning" in caplog.text


def test_default_log(caplog, capsys):
    """Default message only in STDOUT."""
    caplog.set_level(logging.INFO)
    logger = logging.getLogger()

    with setup_logging(logger, level=logging.INFO):
        logger.info("foo")
        logger.debug("A debug message")

    stdout_log = capsys.readouterr().out

    assert "Logging to file is disabled." not in caplog.text  # DEBUG message
    assert "INFO" not in stdout_log
    assert "foo" in stdout_log
    assert "A debug message" not in stdout_log


def test_info_log(caplog, monkeypatch, tmp_path, capsys):
    """Default message in STDOUT and file."""
    monkeypatch.chdir(tmp_path)
    caplog.set_level(logging.INFO)
    logger = logging.getLogger()

    with setup_logging(logger, logfile="logfile.log", level=logging.INFO):
        logger.info("foo")
        logger.debug("A debug message")

    with Path.open("logfile.log") as f:
        file_log = f.read()

    stdout_log = capsys.readouterr().out
    log_path = str((tmp_path / "logfile.log").absolute().resolve())

    assert file_log == stdout_log
    assert f"This log is also available at '{log_path}'" in caplog.text

    for logtext in [stdout_log, file_log]:
        assert "INFO" not in logtext
        assert "foo" in logtext
        assert "A debug message" not in logtext


def test_debug_log(caplog, monkeypatch, tmp_path, capsys):
    """Debug message in STDOUT and file."""
    monkeypatch.chdir(tmp_path)
    caplog.set_level(logging.DEBUG)
    logger = logging.getLogger()

    with setup_logging(logger, logfile="logfile.log", level=logging.DEBUG):
        logger.info("foo")
        logger.debug("A debug message")

    with Path.open("logfile.log") as f:
        file_log = f.read()

    stdout_log = capsys.readouterr().out
    log_path = str((tmp_path / "logfile.log").absolute())

    assert file_log == stdout_log
    assert f"This log is also available at '{log_path}'" in caplog.text

    for logtext in [stdout_log, file_log]:
        assert "foo" in logtext
        assert "A debug message" in logtext
        # Test that debug information is in output
        assert "test_logger.py:test_debug_log:" in logtext
