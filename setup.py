#!/usr/bin/env python3
# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
#
# Copyright (c) 2021 Authors and contributors
#
# Released under the GNU Public Licence, v2 or any higher version
# SPDX-License-Identifier: GPL-2.0-or-later
"""Setuptools-based setup script for tests of mdacli.

For a basic installation just type the command::

  python setup.py develop

"""

import re

# Always prefer setuptools over distutils
from setuptools import setup


VERSION = re.search(
    r'__version__\s*=\s*[\'"]([^\'"]*)[\'"]', open("src/mdacli/__init__.py").read()
    ).group(1)

setup(version=VERSION)
