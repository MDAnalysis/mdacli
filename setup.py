#!/usr/bin/env python3
# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
#
# Copyright (c) 2021 Authors and contributors
#
# Released under the GNU Public Licence, v2 or any higher version
# SPDX-License-Identifier: GPL-2.0-or-later

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
    version='0.0.0',
    description='A command line client for MDAnalysis Analysis classes.',
    long_description=long_description,
    long_description_content_type='text/x-rst',
    license='GPLv2',
    author='PicoCentauri, joaomcteixeira',
    author_email='philip-loche@gmx.de, joaomcteixeira@gmail.com',
    maintainer='PicoCentauri, joaomcteixeira, MDA devs',
    maintainer_email='philip-loche@gmx.de, joaomcteixeira@gmail.com, mdanalysis@numfocus.org',
    url='https://www.mdanalysis.org/',
    packages=find_packages('src'),
    package_dir={'': 'src'},
    py_modules=[splitext(basename(i))[0] for i in glob("src/*.py")],
    include_package_data=True,
    zip_safe=False,
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Environment :: Console',
        'License :: OSI Approved :: GNU General Public License v2 (GPLv2)',
        'Intended Audience :: Science/Research',
        'Natural Language :: English',
        'Operating System :: POSIX',
        'Operating System :: MacOS',
        'Operating System :: Microsoft',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Topic :: Scientific/Engineering :: Chemistry',
        ],
    project_urls={
        'webpage': '',
        'Documentation': '',
        'Changelog': '',
        'Issue Tracker': '',
        'Discussion Forum': '',
        },
    keywords=[
        # eg: 'keyword1', 'keyword2', 'keyword3',
        ],
    python_requires='>=3.6, <3.10',
    install_requires=[
        # 'click',
        # eg: 'aspectlib==1.1.1', 'six>=1.7',
        ],
    extras_require={
        # eg:
        #   'rst': ['docutils>=0.11'],
        #   ':python_version=="2.6"': ['argparse'],
        },
    setup_requires=[
        #   'pytest-runner',
        #   'setuptools_scm>=3.3.1',
        ],
    entry_points={
        'console_scripts': [
            'mdacli= mdacli:main',
            ]
        #
        },
    # cmdclass={'build_ext': optional_build_ext},
    # ext_modules=[
    #    Extension(
    #        splitext(relpath(path, 'src').replace(os.sep, '.'))[0],
    #        sources=[path],
    #        include_dirs=[dirname(path)]
    #    )
    #    for root, _, _ in os.walk('src')
    #    for path in glob(join(root, '*.c'))
    # ],
    )
