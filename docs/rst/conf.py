# -*- coding: utf-8 -*-
"""Config file for Sphinx-docs."""
from __future__ import unicode_literals

import os
import mock
import sys

import msmb_theme
import sphinx_rtd_theme


# activate if there are dependencies
mock_modules = [
    'matplotlib',
    'MDAnalysis',
    'MDAnalysis.analysis',
    'MDAnalysis.analysis.base',
    'MDAnalysis.transformations',
    'MDAnalysis.transformations.boxdimensions',
    'numpy',
    ]

for modulename in mock_modules:
    sys.modules[modulename] = mock.Mock()

extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.autosummary',
    'sphinx.ext.coverage',
    'sphinx.ext.doctest',
    'sphinx.ext.extlinks',
    'sphinx.ext.ifconfig',
    'sphinx.ext.intersphinx',
    'sphinx.ext.napoleon',
    'sphinx.ext.todo',
    'sphinx.ext.viewcode',
    'sphinx.ext.autosectionlabel',
    'sphinx_rtd_theme',
    'sphinx_sitemap',
    ]

# for sitemap with https://github.com/jdillard/sphinx-sitemap
site_url = "https://mdacli.mdanalysis.org/"

todo_include_todos = True

exclude_patterns = [
    'i_*',
    ]

if os.getenv('SPELLCHECK'):
    extensions += 'sphinxcontrib.spelling',
    spelling_show_suggestions = True
    spelling_lang = 'en_US'
    # https://sphinxcontrib-spelling.readthedocs.io/en/latest/customize.html
    spelling_word_list_filename = ['../spelling_wordlist.txt']

source_suffix = '.rst'
master_doc = 'index'
project = 'mdacli'
year = '2021'
author = 'Philip Loche, Joao MC Teixeira and Oliver Beckstein'
copyright = '{0}, {1}'.format(year, author)
version = release = '0.1.8'

pygments_style = 'trac'
templates_path = ['_templates']

extlinks = {
    'issue': ('https://github.com/MDAnalysis/mdacli/cissues/%s', '#'),
    'pr': ('https://github.com/MDAnalysis/mdacli/pull/%s', 'PR #'),
    }

# codecov io closes connection if host is accessed too repetitively.
# codecov links are ignored here for the same reason there's a sleep
# in the .travis.yml file
# see https://github.com/codecov/codecov-python/issues/158
linkcheck_ignore = [
    r'https://codecov.io/gh/MDAnalysis/mdacli/*',
    ]

html_theme = "msmb_theme"
html_theme_path = [msmb_theme.get_html_theme_path(),
                   sphinx_rtd_theme.get_html_theme_path()]
html_theme_options = {
    'canonical_url': '',
    'logo_only': True,
    'display_version': True,
    'prev_next_buttons_location': 'bottom',
    'style_external_links': False,
    'style_nav_header_background': 'white',
    # Toc options
    'collapse_navigation': True,
    'sticky_navigation': True,
    'navigation_depth': 4,
    'includehidden': True,
    'titles_only': False,
}
html_use_smartypants = True
html_last_updated_fmt = '%b %d, %Y'
html_split_index = False
html_short_title = '%s-%s' % (project, version)
html_baseurl = site_url
html_logo = "_static/logos/mdacli-logo.png"
html_favicon = "_static/logos/mdanalysis-logo.ico"
html_static_path = ['_static']
html_css_files = ['custom.css']
html_use_opensearch = site_url

html_context = {
    'versions_json_url': os.path.join(site_url, "versions.json")
}

napoleon_use_ivar = True
napoleon_use_rtype = False
napoleon_use_param = False

# Configuration for intersphinx: refer to the Python standard library
# and other packages used by mdacli
intersphinx_mapping = {'https://docs.python.org/': None,
                       'https://docs.mdanalysis.org/stable/': None,
                       }
