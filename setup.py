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

from glob import glob
from os.path import basename, dirname, join, splitext

from setuptools import find_packages, setup


def read(*names, **kwargs):
    """Read description files."""
    path = join(dirname(__file__), *names)
    with open(path, encoding=kwargs.get('encoding', 'utf8')) as fh:
        return fh.read()


long_description = '{}\n{}'.format(
    read('README.rst'),
    read(join('docs', 'CHANGELOG.rst')),
    )

setup(
    name='mdacli',
    version='0.1.0',
    description='A command line client for MDAnalysis Analysis classes.',
    long_description=long_description,
    long_description_content_type='text/x-rst',
    license='GPLv3',
    author='PicoCentauri, joaomcteixeira',
    author_email='philip-loche@gmx.de, joaomcteixeira@gmail.com',
    maintainer='PicoCentauri, joaomcteixeira, MDA devs',
    maintainer_email='philip-loche@gmx.de, joaomcteixeira@gmail.com, mdanalysis@numfocus.org',
    url='https://github.com/MDAnalysis/mdacli',
    packages=find_packages('src'),
    package_dir={'': 'src'},
    py_modules=[splitext(basename(i))[0] for i in glob("src/*.py")],
    include_package_data=True,
    zip_safe=False,
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Environment :: Console',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        'Intended Audience :: Science/Research',
        'Natural Language :: English',
        'Operating System :: POSIX',
        'Operating System :: MacOS :: MacOS X',
        'Operating System :: Microsoft :: Windows ',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Topic :: Scientific/Engineering',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Topic :: Scientific/Engineering :: Chemistry',
        'Topic :: Scientific/Engineering :: Physics',
        'Topic :: Software Development :: Libraries :: Python Modules',
        'Topic :: System :: Shells',
        ],
    project_urls={
        'webpage': '',
        'Documentation': '',
        'Changelog': '',
        'Issue Tracker': '',
        'Discussion Forum': '',
        },
    keywords=[
        'Science',
        'Molecular Dynamics',
        'MDAnalysis',
    ],
    python_requires='>=3.7',
    install_requires=[
            'MDAnalysis>=2.0.0',
        ],
    entry_points={
        'console_scripts': [
            'mda = mdacli.cli:main',
            ]
        },
