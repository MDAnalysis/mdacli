MDAnalysis command line interface
=================================

|pypi| |mdanalysis| |codecov| |docs| |test|

``mdacli`` is a simple command line interface (CLI) to the analysis classes of `MDAnalysis`_
using argparse_. This project is in an **early development stage** and
work in progress. `Contributions are welcome <https://github.com/MDAnalysis/mdacli/blob/main/docs/CONTRIBUTING.rst>`_!

To install `mdacli` refer to the `INSTALL file <https://github.com/MDAnalysis/mdacli/blob/main/docs/rst/installation.rst>`_.

Run `mdacli`::

   mda -h

For a help and an overview of the supported modules. A help
message for each module is available using::

   mda <module> -h


Available modules
-----------------

Currently the following analysis modules are available

.. list-table::
   :widths: 25 50
   :header-rows: 1

   * - Module Name
     - Description

   * - AlignTraj
     - RMS-align trajectory to a reference structure using a selection.
   * - AverageStructure
     - RMS-align trajectory to a reference structure using a selection,
       and calculate the average coordinates of the trajectory.
   * - Contacts
     - Calculate contacts based observables.
   * - DensityAnalysis
     - Volumetric density analysis.
   * - DistanceMatrix
     - Calculate the pairwise distance between each frame in a trajectory
   * - Dihedral
     - Calculate dihedral angles for specified atomgroups.
   * - Janin
     - Calculate χ_1 and χ_2 dihedral angles of selected group
   * - Ramachandran
     - Calculate ϕ and ψ dihedral angles of selected group
   * - DielectricConstant
     - Computes the average dipole moment.
   * - GNMAnalysis
     - Basic tool for GNM analysis.
   * - closeContactGNMAnalysis
     - GNMAnalysis only using close contacts.
   * - HELANAL
     - Perform HELANAL helix analysis on your trajectory.
   * - HoleAnalysis
     - Run `hole` program on a trajectory.
   * - LinearDensity
     - Linear density profile
   * - EinsteinMSD
     - Class to calculate Mean Squared Displacement by the Einstein relation.
   * - PCA
     - Principal component analysis on an MD trajectory.
   * - InterRDF
     - Intermolecular pair distribution function
   * - RMSD
     - Class to perform RMSD analysis on a trajectory.
   * - RMSF
     - Calculate RMSF of given atoms across a trajectory.

More information about each module is available through the help
page or at the `MDAnalysis documentation`_.

.. _argparse: https://docs.python.org/3/library/argparse.html
.. _MDAnalysis: https://www.mdanalysis.org
.. _`MDAnalysis installed`: https://userguide.mdanalysis.org/stable/installation.html
.. _`MDAnalysis documentation`: https://docs.mdanalysis.org/stable/documentation_pages/analysis_modules.html

 .. |pypi| image:: https://img.shields.io/pypi/v/mdacli.svg
   :alt: PyPI Package latest release
   :target: https://pypi.org/project/mdacli

 .. |mdanalysis| image:: https://img.shields.io/badge/powered%20by-MDAnalysis-orange.svg?logoWidth=16&logo=data:image/x-icon;base64,AAABAAEAEBAAAAEAIAAoBAAAFgAAACgAAAAQAAAAIAAAAAEAIAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAJD+XwCY/fEAkf3uAJf97wGT/a+HfHaoiIWE7n9/f+6Hh4fvgICAjwAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACT/yYAlP//AJ///wCg//8JjvOchXly1oaGhv+Ghob/j4+P/39/f3IAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAJH8aQCY/8wAkv2kfY+elJ6al/yVlZX7iIiI8H9/f7h/f38UAAAAAAAAAAAAAAAAAAAAAAAAAAB/f38egYF/noqAebF8gYaagnx3oFpUUtZpaWr/WFhY8zo6OmT///8BAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAgICAn46Ojv+Hh4b/jouJ/4iGhfcAAADnAAAA/wAAAP8AAADIAAAAAwCj/zIAnf2VAJD/PAAAAAAAAAAAAAAAAICAgNGHh4f/gICA/4SEhP+Xl5f/AwMD/wAAAP8AAAD/AAAA/wAAAB8Aov9/ALr//wCS/Z0AAAAAAAAAAAAAAACBgYGOjo6O/4mJif+Pj4//iYmJ/wAAAOAAAAD+AAAA/wAAAP8AAABhAP7+FgCi/38Axf4fAAAAAAAAAAAAAAAAiIiID4GBgYKCgoKogoB+fYSEgZhgYGDZXl5e/m9vb/9ISEjpEBAQxw8AAFQAAAAAAAAANQAAADcAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAjo6Mb5iYmP+cnJz/jY2N95CQkO4pKSn/AAAA7gAAAP0AAAD7AAAAhgAAAAEAAAAAAAAAAACL/gsAkv2uAJX/QQAAAAB9fX3egoKC/4CAgP+NjY3/c3Nz+wAAAP8AAAD/AAAA/wAAAPUAAAAcAAAAAAAAAAAAnP4NAJL9rgCR/0YAAAAAfX19w4ODg/98fHz/i4uL/4qKivwAAAD/AAAA/wAAAP8AAAD1AAAAGwAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAALGxsVyqqqr/mpqa/6mpqf9KSUn/AAAA5QAAAPkAAAD5AAAAhQAAAAEAAAAAAAAAAAAAAAAAAAAAAAAAAAAAADkUFBSuZ2dn/3V1df8uLi7bAAAATgBGfyQAAAA2AAAAMwAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAB0AAADoAAAA/wAAAP8AAAD/AAAAWgC3/2AAnv3eAJ/+dgAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA9AAAA/wAAAP8AAAD/AAAA/wAKDzEAnP3WAKn//wCS/OgAf/8MAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAIQAAANwAAADtAAAA7QAAAMAAABUMAJn9gwCe/e0Aj/2LAP//AQAAAAAAAAAA
   :alt: Powered by MDAnalysis
   :target: https://www.mdanalysis.org

 .. |docs| image:: https://readthedocs.org/projects/mdacli/badge/?version=latest
   :target: https://mdacli.readthedocs.io/en/latest/?badge=latest
   :alt: Documentation Status

 .. |test| image:: https://github.com/MDAnalysis/mdacli/actions/workflows/test.yml/badge.svg?branch=main
   :alt: Github Actions Test Status
   :target: https://github.com/MDAnalysis/mdacli/actions/workflows/test.yml

 .. |codecov| image:: https://codecov.io/gh/MDAnalysis/mdacli/branch/main/graph/badge.svg?token=ets2mZ6xJD
    :alt: Codecov mdacli
    :target: https://codecov.io/gh/MDAnalysis/mdacli
