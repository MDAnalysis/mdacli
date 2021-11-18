============
Installation
============

``mdacli`` is a *pluging* tool for `MDAnalysis
<https://www.mdanalysis.org/>`_.  So the installation of `mdacli`
depends slightly on your installation of `MDAnalysis`.  If you have
`MDAnalysis v2` already installed, you should install `mdacli` on top of
it in the same Python environment. So, activate that environment first.

If you just want to use `mdacli`, you can install it from `PyPI`::

    pip install mdacli --no-deps

Or if you are a `github` expert you can clone the repository and install
it from the repository directly::

    python setup.py develop --no-deps

If you don't have `MDAnalysis v2` installed and want to install
everything together::

    pip install mdacli

This should install also MDAnalysis version 2.

Any doubts please `contact us <https://github.com/MDAnalysis/mdacli/issues>`.
