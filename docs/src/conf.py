"""Config file for Sphinx-docs."""

import os
import sys
from unittest import mock

# activate if there are dependencies
mock_modules = [
    "matplotlib",
    "MDAnalysis",
    "MDAnalysis.analysis",
    "MDAnalysis.analysis.base",
    "MDAnalysis.transformations",
    "MDAnalysis.transformations.boxdimensions",
    "numpy",
]

for modulename in mock_modules:
    sys.modules[modulename] = mock.Mock()

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    "sphinx.ext.coverage",
    "sphinx.ext.doctest",
    "sphinx.ext.extlinks",
    "sphinx.ext.ifconfig",
    "sphinx.ext.intersphinx",
    "sphinx.ext.napoleon",
    "sphinx.ext.todo",
    "sphinx.ext.viewcode",
    "sphinx.ext.autosectionlabel",
    "mdanalysis_sphinx_theme",
    "sphinx_sitemap",
]

# for sitemap with https://github.com/jdillard/sphinx-sitemap
site_url = "https://mdacli.mdanalysis.org/"

todo_include_todos = True

exclude_patterns = [
    "i_*",
]

if os.getenv("SPELLCHECK"):
    extensions += ("sphinxcontrib.spelling",)
    spelling_show_suggestions = True
    spelling_lang = "en_US"
    # https://sphinxcontrib-spelling.readthedocs.io/en/latest/customize.html
    spelling_word_list_filename = ["../spelling_wordlist.txt"]

source_suffix = ".rst"
master_doc = "index"
project = "mdacli"
year = "2021"
author = "Philip Loche, Joao MC Teixeira and Oliver Beckstein"
copyright = f"{year}, {author}"
version = release = "0.1.32"

pygments_style = "trac"
templates_path = ["_templates"]

extlinks = {
    "issue": ("https://github.com/MDAnalysis/mdacli/cissues/%s", "#"),
    "pr": ("https://github.com/MDAnalysis/mdacli/pull/%s", "PR #"),
}

# codecov io closes connection if host is accessed too repetitively.
# codecov links are ignored here for the same reason there's a sleep
# in the .travis.yml file
# see https://github.com/codecov/codecov-python/issues/158
linkcheck_ignore = [
    r"https://codecov.io/gh/MDAnalysis/mdacli/*",
]

html_theme = "mdanalysis_sphinx_theme"
htmlhelp_basename = "mdacli"

html_theme_options = {
    "mda_official": False,
}
html_use_smartypants = True
html_last_updated_fmt = "%b %d, %Y"
html_split_index = False
html_short_title = f"{project}-{version}"
html_baseurl = site_url
html_logo = "_static/logos/mdacli-logo.png"
html_favicon = "_static/logos/mdanalysis-logo.ico"
html_static_path = ["_static"]
html_css_files = ["custom.css"]
html_use_opensearch = site_url

html_context = {"versions_json_url": f"{site_url}/versions.json"}

napoleon_use_ivar = True
napoleon_use_rtype = False
napoleon_use_param = False

# Configuration for intersphinx: refer to the Python standard library
# and other packages used by mdacli
intersphinx_mapping = {
    "python": ("https://docs.python.org/3/", None),
    "MDAnalysis": ("https://docs.mdanalysis.org/stable/", None),
}
