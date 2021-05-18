#!/usr/bin/env python3
# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
#
# Copyright (c) 2021 Authors and contributors
#
# Released under GNU Public Licence, v2 or any higher version
# SPDX-License-Identifier: GPL-2.0-or-later

"""
Emphaising strings with colors etc.

Taken from https://gist.github.com/tuvokki/14deb97bef6df9bc6553.
"""


class Emphasise:
    """Class for emphaising strings with colors etc.

    Attributes
    ----------
    bold : str
        bold attribute
    underline : str
        underline attribute
    gray : str
        gray color
    red : str
        red color
    green : str
        green color
    yellow : str
        yellow color
    blue : str
        blue color
    pink : str
        pink color
    turquoise : str
        turquoise color
    """

    _endc = '\033[0m'
    bold = '\033[1m'
    underline = '\033[4m'
    gray = '\033[90m'
    red = '\033[91m'
    green = '\033[92m'
    yellow = '\033[93m'
    blue = '\033[94m'
    pink = '\033[95m'
    turquoise = '\033[96m'

    def _emphasise(self, str, style):
        return f"{style}{str}{self._endc}"

    @staticmethod
    def emphasise(str, style):
        """Decorate a ``str`` with desired style.

        The Style could be a color, bold or underline.

        Parameters
        ----------
        message : str
            message to print
        style : str
            emphasising style. See class attributes for available styles

        Example
        -------
        >>> print(Emphasise.emphasise("My colored message", Emphasise.blue))
        """
        style_str = getattr(Emphasise, style)
        return Emphasise._emphasise(Emphasise, str, style_str)

    @staticmethod
    def warning(message):
        """Return a yellow warning.

        Parameters
        ----------
        message : str
            yellow warning to print

        Returns
        -------
        decorated_message : str
            decorated message

        Example
        ------
        >>> print(bcolors.warning("Potential Danger!"))
        """
        return Emphasise._emphasise(Emphasise, message, Emphasise.yellow)

    @staticmethod
    def error(message):
        """Return a red error.

        Parameters
        ----------
        message : str
            red error to return

        Returns
        -------
        decorated_message : str
            decorated message

        Example
        ------
        >>> print(bcolors.error("Potential Danger!"))
        """
        return Emphasise._emphasise(Emphasise, message, Emphasise.red)

    @staticmethod
    def ok(message):
        """Return a green ok.

        Parameters
        ----------
        message : str
            green ok to return

        Returns
        -------
        decorated_message : str
            decorated message

        Example
        ------
        >>> print(bcolors.ok("Yay!"))
        """
        return Emphasise._emphasise(Emphasise, message, Emphasise.green)

    @staticmethod
    def info(message):
        """Return a blue info.

        Parameters
        ----------
        message : str
            blue info to return

        Returns
        -------
        decorated_message : str
            decorated message

        Example
        ------
        >>> print(bcolors.info("Blue Yay!"))
        """
        return Emphasise._emphasise(Emphasise, message, Emphasise.blue)

    @staticmethod
    def header(message):
        """Return a pink header.

        Parameters
        ----------
        message : str
            pink header to return

        Returns
        -------
        decorated_message : str
            decorated message

        Example
        ------
        >>> print(bcolors.header("This is great"))
        """
        return Emphasise._emphasise(Emphasise, message, Emphasise.blue)

    @staticmethod
    def debug(message):
        """Return a turquoise debug message.

        Parameters
        ----------
        message : str
            turquoise debug message to return

        Returns
        -------
        decorated_message : str
            decorated message

        Example
        ------
        >>> print(bcolors.debug("a=1"))
        """
        return Emphasise._emphasise(Emphasise, message, Emphasise.turquoise)
