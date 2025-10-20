#!/usr/bin/env python3
#
# Copyright (c) 2021 Authors and contributors
#
# Released under the GNU Public Licence, v2 or any higher version
# SPDX-License-Identifier: GPL-2.0-or-later
"""Logging."""

import contextlib
import logging
import sys
import warnings
from pathlib import Path

from .colors import Emphasise


def check_suffix(filename: str | Path, suffix: str) -> str | Path:
    """Check the suffix of a file name and adds if it not existing.

    If ``filename`` does not end with ``suffix`` the ``suffix`` is added and a
    warning will be issued.

    Parameters
    ----------
    filename : Name of the file to be checked.
    suffix : Expected file suffix i.e. ``.txt``.

    Returns
    -------
    Checked and probably extended file name.
    """
    path_filename = Path(filename)

    if path_filename.suffix != suffix:
        warnings.warn(
            f"The file name should have a '{suffix}' extension. The user "
            f"requested the file with name '{filename}', but it will be saved "
            f"as '{filename}{suffix}'.",
            stacklevel=1,
        )
        path_filename = path_filename.parent / (path_filename.name + suffix)

    if type(filename) is str:
        return str(path_filename)
    return path_filename


@contextlib.contextmanager
def setup_logging(
    logobj: logging.Logger,
    logfile: str | Path | None = None,
    level: int = logging.WARNING,
):
    """Create a logging environment for a given ``log_obj``.

    Parameters
    ----------
    logobj : ``logging.Logger``
        A logging instance
    logfile : str
        Name of the log file
    level : int
        Set the root logger level to the specified level. If for example set
        to :py:obj:`logging.DEBUG` detailed debug logs inludcing filename and
        function name are displayed. For :py:obj:`logging.INFO` only the
        message logged from errors, warnings and infos will be displayed.
    """
    try:
        format = ""
        if level == logging.DEBUG:
            format += "[{levelname}]:{filename}:{funcName}:{lineno} - "
        format += "{message}"

        formatter = logging.Formatter(format, style="{")
        handlers: list[logging.StreamHandler | logging.FileHandler] = []

        stream_handler = logging.StreamHandler(sys.stdout)
        stream_handler.setFormatter(formatter)
        handlers.append(stream_handler)

        if logfile:
            logfile = check_suffix(filename=logfile, suffix=".log")
            file_handler = logging.FileHandler(filename=str(logfile), encoding="utf-8")
            file_handler.setFormatter(formatter)
            handlers.append(file_handler)
        else:
            logging.addLevelName(logging.INFO, Emphasise.info("INFO"))
            logging.addLevelName(logging.DEBUG, Emphasise.debug("DEBUG"))
            logging.addLevelName(logging.WARNING, Emphasise.warning("WARNING"))
            logging.addLevelName(logging.ERROR, Emphasise.error("ERROR"))

        logging.basicConfig(format=format, handlers=handlers, level=level, style="{")
        logging.captureWarnings(True)

        if logfile:
            abs_path = str(Path(logfile).absolute().resolve())
            logobj.info(f"This log is also available at '{abs_path}'.")
        else:
            logobj.debug("Logging to file is disabled.")

        for handler in handlers:
            logobj.addHandler(handler)

        yield

    finally:
        for handler in handlers:
            handler.flush()
            handler.close()
            logobj.removeHandler(handler)
