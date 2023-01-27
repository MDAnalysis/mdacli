
Changelog
=========


v0.1.24 (2023-01-27)
------------------------------------------

* Fix classifier in setup.cfg

v0.1.23 (2023-01-27)
------------------------------------------

* Disable docs deploy on Github Pages
* Add Python 3.11 to CI matrix

v0.1.22 (2023-01-27)
------------------------------------------

* Limit number of threads used by analysis

v0.1.21 (2022-04-22)
------------------------------------------

* Allow None as start/stop values (need for 1 frame trajectories)

v0.1.20 (2022-04-13)
------------------------------------------

* Fix isort dependency
* Added tests for Python 3.10
* Introduce new setup.cfg and pyproject.toml
* Update __authors__

v0.1.19 (2022-04-12)
------------------------------------------

* Update README.rst for MDAnalysis 2.1.0 modules

v0.1.18 (2022-04-12)
------------------------------------------

* Fixed typo in libcli.py

v0.1.17 (2022-04-07)
------------------------------------------

* Added dihedral module

v0.1.16 (2022-02-25)
------------------------------------------

* Do not convert None types

v0.1.15 (2022-02-25)
------------------------------------------

* Set positional arguments as required in cli

v0.1.14 (2022-02-17)
------------------------------------------

* corrects bump2version changelog update

v0.1.13 (2022-02-16)
------------------------------------------

* Added conda package install instructions (#88)

v0.1.12 (2022-01-19)
-------------------------------------------------------------------------

* Support list of AtomGroups as parameters (#82)
* Simplify `add_argument` logic in `create_CLI` (#82)
* Allow list of reference classes in module detection (#82)
* Support for generic classes as reference in module detection (#82)
* Rename `save_results`` to `save` (#82)
* More tests for docstring parsing and CLI creation (#82)

v0.1.11 (2022-01-19)
-------------------------------------------------------------------------

* Improved help for run parameters (#83)

v0.1.10 (2022-01-18)
------------------------------------------

* Removed conda dependency from CI and tox (#86)

v0.1.9 (2022-01-16)
------------------------------------------

* Fix test banner in README.rst (#85)

v0.1.8 (2022-01-16)
------------------------------------------

* Use Github actions matrix for tests (#68)
* Fix Conda permissions on MacOS (#68)
* Fix Tests failing on Windows (#68)

v0.1.7 (2021-12-18)
------------------------------------------

* Improves regex to convert from time to frame (#81)

v0.1.6 (2021-12-01)
-------------------------------------------

* Fixed URL in docs (#80)

v0.1.5 (2021-12-01)
--------------------------------------------------

* Add doc deployment to CI (#78)

v0.1.4 (2021-11-24)
-------------------------------------------------------------------------

* Link docs to mdacli.mdanalysis.org (#75)

v0.1.3 (2021-11-24)
------------------------------------------

* MDA-style documentation pages (#70)

v0.1.2 (2021-11-18)
------------------------------------------

* Added option to manually set box dimensions (#65)

v0.1.1 (2021-11-18)
------------------------------------------

* corrects .bumpversion.cfg for CHANGELOG
* updates docs/CONTRIBUTING.rst accordingly

v0.1.0 (2021-11-18)
-------------------
* Initial release
