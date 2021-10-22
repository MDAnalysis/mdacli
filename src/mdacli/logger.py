#!/usr/bin/env python3
# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
#
# Copyright (c) 2021 Authors and contributors
#
# Released under the GNU Public Licence, v2 or any higher version
# SPDX-License-Identifier: GPL-2.0-or-later
"""Logging."""

import contextlib
import logging
import sys

from .colors import Emphasise


@contextlib.contextmanager
def setup_logging(logobj, logfile=None, debug=False):
    """
    Create a logging environment for a given logobj.

    Parameters
    ----------
    logobj : ``logging.Logger``
        A logging instance
    logfile : str
        Name of the log file
    debug : bool
        If ``True`` detailed debug logs inludcing filename and function name
        are displayed. If ``False`` only the message logged from
        errors, warnings and infos will be displayed.
    """
    try:
        format = '{message}'
        if debug:
            format = "[{levelname}] {filename}:{name}:{funcName}:{lineno}: " \
                     + format
            level = logging.DEBUG
        else:
            level = logging.INFO

        logging.basicConfig(format=format,
                            handlers=[logging.StreamHandler(sys.stdout)],
                            level=level,
                            style='{')

        if logfile:
            logfile += ".log" * (not logfile.endswith("log"))
            handler = logging.FileHandler(filename=logfile, encoding='utf-8')
            handler.setFormatter(logging.Formatter(format, style='{'))
            logobj.addHandler(handler)
        else:
            logging.addLevelName(logging.INFO, Emphasise.info("INFO"))
            logging.addLevelName(logging.DEBUG, Emphasise.debug("DEBUG"))
            logging.addLevelName(logging.WARNING, Emphasise.warning("WARNING"))
            logging.addLevelName(logging.ERROR, Emphasise.error("ERROR"))
            logobj.info('Logging to file is disabled.')

        yield
    finally:
        handlers = logobj.handlers[:]
        for handler in handlers:
            handler.close()
            logobj.removeHandler(handler)
