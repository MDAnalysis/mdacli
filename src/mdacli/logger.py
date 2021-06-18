#!/usr/bin/env python3
# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
#
# Copyright (c) 2021 Authors and contributors
#
# Released under the GNU Public Licence, v2 or any higher version
# SPDX-License-Identifier: GPL-2.0-or-later

import contextlib
import logging

from logging.handlers import RotatingFileHandler

from mdacli.colors import Emphasise

DEBUGFORMATTER = '{levelname}:{filename}:{name}:{funcName}:%{lineno}: {message}'
"""Debug file formatter."""

INFOFORMATTER = '{levelname}:{message}'
"""Log file and stream output formatter."""

logger = logging.getLogger(__name__)

@contextlib.contextmanager
def setup_logging(logfile=None, debug=False):
    """Setup logging."""
    try:
        log = logging.getLogger()
        if debug:
            logging.basicConfig(format=DEBUGFORMATTER, 
                                level=logging.DEBUG,
                                style='{')
            fmt = logging.Formatter(DEBUGFORMATTER, style='{')
        else:
            logging.basicConfig(format=INFOFORMATTER,
                                level=logging.INFO,
                                style='{')
            fmt = logging.Formatter(INFOFORMATTER, style='{')

        if logfile:
            logfile += ".log" * (not logfile.endswith("log"))
            handler = RotatingFileHandler(filename=logfile,
                                          encoding='utf-8')
            handler.setFormatter(fmt)
            log.addHandler(handler)
        else:
            logging.addLevelName(logging.INFO, Emphasise.info("INFO"))
            logging.addLevelName(logging.DEBUG, Emphasise.debug("DEBUG"))
            logging.addLevelName(logging.WARNING, Emphasise.warning("ERROR"))
            logging.addLevelName(logging.ERROR, Emphasise.error("ERROR"))
            logger.info('Logging to file is disabled')

        yield
    finally:
        handlers = log.handlers[:]
        for handler in handlers:
            handler.close()
            log.removeHandler(handler)
