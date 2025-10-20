============
Installation
============

``mdacli`` is a *pluging* tool for `MDAnalysis
<https://www.mdanalysis.org/>`_.  So the installation of `mdacli`
depends slightly on your installation of `MDAnalysis`.  If you have
`MDAnalysis` already installed, you should install `mdacli` on top of
it in the same Python environment. So, activate that environment first.

pip
---

If you just want to use `mdacli`, you can install it with `pip`_::

    pip install mdacli --no-deps --upgrade

Or if you are a `Github` expert you can clone the repository and install
it from the repository directly::

    python setup.py develop --no-deps

If you don't have `MDAnalysis` installed and want to install
everything together::

    pip install mdacli --upgrade

conda
-----

First installation with conda_::

    conda install -c conda-forge mdacli

which will automatically also install `MDAnalysis`.

To upgrade later::

   conda update mdacli

Any doubts please `contact us <https://github.com/MDAnalysis/mdacli/issues>`.

.. _pip:
   http://www.pip-installer.org/en/latest/index.html
.. _conda:
   http://conda.pydata.org/docs/
