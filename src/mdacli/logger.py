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

from mdacli.colors import Emphasise


DEBUGFORMATTER = '[{levelname}] {filename}:{name}:{funcName}:{lineno}:'
DEBUGFORMATTER += '{message}'
"""Debug file formatter."""

INFOFORMATTER = '{message}'
"""Log file and stream output formatter."""

logger = logging.getLogger(__name__)


@contextlib.contextmanager
def setup_logging(logfile=None, debug=False):
    """
    Create a logging environment.

    Parameters
    ----------
    logfile : str
        Name of the log file
    debug : bool
        Display debug logs. If ``False`` error, warnings and infos will be
        logged.
    """
    try:
        log = logging.getLogger()
        if debug:
            format = DEBUGFORMATTER
            level = logging.DEBUG
        else:
            format = INFOFORMATTER
            level = logging.INFO

        logging.basicConfig(format=format,
                            handlers=[logging.StreamHandler(sys.stdout)],
                            level=level,
                            style='{')

        if logfile:
            logfile += ".log" * (not logfile.endswith("log"))
            handler = logging.FileHandler(filename=logfile, encoding='utf-8')
            handler.setFormatter(logging.Formatter(format, style='{'))
            log.addHandler(handler)
        else:
            logging.addLevelName(logging.INFO, Emphasise.info("INFO"))
            logging.addLevelName(logging.DEBUG, Emphasise.debug("DEBUG"))
            logging.addLevelName(logging.WARNING, Emphasise.warning("WARNING"))
            logging.addLevelName(logging.ERROR, Emphasise.error("ERROR"))
            logger.info('Logging to file is disabled')

        yield
    finally:
        handlers = log.handlers[:]
        for handler in handlers:
            handler.close()
            log.removeHandler(handler)
