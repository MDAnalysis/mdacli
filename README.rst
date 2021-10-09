MDAnalysis command line interface README
========================================

|mdanalysis| 

|test_py36| |test_py37| |test_py38| |test_py39| 

A simple command line interface (CLI) to the analysis classes of `MDAnalysis`_
using argparse_. This project is in an **early development stage** and 
work in progress. Contributions are welcome!

Provided you have `MDAnalysis installed`_, get the mdacli source from 
https://github.com/MDAnalysis/mdacli and install ``mdacli`` with the following::

   python setup.py develop

and run::

   mda -h

for a help and an overview of the supported modules. A help 
message for each module is available using::

   mda module -h


Currentlty the following analysis modules are available

.. list-table::
   :widths: 25 50
   :header-rows: 1

   * - Module Name
     - Description
  
   * - AlignTraj           
     - RMS-align trajectory to a reference structure using a selection.
   * - AverageStructure    
     - RMS-align trajectory to a reference structure using a selection,
   * - Contacts            
     - Calculate contacts based observables.
   * - DensityAnalysis     
     - Volumetric density analysis.
   * - DistanceMatrix      
     - Calculate the pairwise distance between each frame in a trajectory
   * - Janin               
     - Calculate :math:`\chi_1` and :math:`\chi_2` dihedral angles of selected 
       group
   * - Ramachandran        
     - Calculate :math:`\phi` and :math:`\psi` dihedral angles of selected group
   * - GNMAnalysis         
     - Basic tool for GNM analysis.
   * - closeContactGNMAnalysis
     - GNMAnalysis only using close contacts.
   * - HELANAL             
     - Perform HELANAL helix analysis on your trajectory.
   * - HoleAnalysis        
     - Run program `hole` on a trajectory.
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

 .. |mdanalysis| image:: https://img.shields.io/badge/powered%20by-MDAnalysis-orange.svg?logoWidth=16&logo=data:image/x-icon;base64,AAABAAEAEBAAAAEAIAAoBAAAFgAAACgAAAAQAAAAIAAAAAEAIAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAJD+XwCY/fEAkf3uAJf97wGT/a+HfHaoiIWE7n9/f+6Hh4fvgICAjwAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACT/yYAlP//AJ///wCg//8JjvOchXly1oaGhv+Ghob/j4+P/39/f3IAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAJH8aQCY/8wAkv2kfY+elJ6al/yVlZX7iIiI8H9/f7h/f38UAAAAAAAAAAAAAAAAAAAAAAAAAAB/f38egYF/noqAebF8gYaagnx3oFpUUtZpaWr/WFhY8zo6OmT///8BAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAgICAn46Ojv+Hh4b/jouJ/4iGhfcAAADnAAAA/wAAAP8AAADIAAAAAwCj/zIAnf2VAJD/PAAAAAAAAAAAAAAAAICAgNGHh4f/gICA/4SEhP+Xl5f/AwMD/wAAAP8AAAD/AAAA/wAAAB8Aov9/ALr//wCS/Z0AAAAAAAAAAAAAAACBgYGOjo6O/4mJif+Pj4//iYmJ/wAAAOAAAAD+AAAA/wAAAP8AAABhAP7+FgCi/38Axf4fAAAAAAAAAAAAAAAAiIiID4GBgYKCgoKogoB+fYSEgZhgYGDZXl5e/m9vb/9ISEjpEBAQxw8AAFQAAAAAAAAANQAAADcAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAjo6Mb5iYmP+cnJz/jY2N95CQkO4pKSn/AAAA7gAAAP0AAAD7AAAAhgAAAAEAAAAAAAAAAACL/gsAkv2uAJX/QQAAAAB9fX3egoKC/4CAgP+NjY3/c3Nz+wAAAP8AAAD/AAAA/wAAAPUAAAAcAAAAAAAAAAAAnP4NAJL9rgCR/0YAAAAAfX19w4ODg/98fHz/i4uL/4qKivwAAAD/AAAA/wAAAP8AAAD1AAAAGwAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAALGxsVyqqqr/mpqa/6mpqf9KSUn/AAAA5QAAAPkAAAD5AAAAhQAAAAEAAAAAAAAAAAAAAAAAAAAAAAAAAAAAADkUFBSuZ2dn/3V1df8uLi7bAAAATgBGfyQAAAA2AAAAMwAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAB0AAADoAAAA/wAAAP8AAAD/AAAAWgC3/2AAnv3eAJ/+dgAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA9AAAA/wAAAP8AAAD/AAAA/wAKDzEAnP3WAKn//wCS/OgAf/8MAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAIQAAANwAAADtAAAA7QAAAMAAABUMAJn9gwCe/e0Aj/2LAP//AQAAAAAAAAAA
    :alt: Powered by MDAnalysis
    :target: https://www.mdanalysis.org

 .. |test_py36| image:: https://github.com/MDAnalysis/mdacli/actions/workflows/py36.yml/badge.svg
   :alt: Github Actions Build Status
   :target: https://github.com/MDAnalysis/mdacli/actions/workflows/py36.ym

 .. |test_py37| image:: https://github.com/MDAnalysis/mdacli/actions/workflows/py37.yml/badge.svg
   :alt: Github Actions Build Status
   :target: https://github.com/MDAnalysis/mdacli/actions/workflows/py37.ym

 .. |test_py38| image:: https://github.com/MDAnalysis/mdacli/actions/workflows/py38.yml/badge.svg
   :alt: Github Actions Build Status
   :target: https://github.com/MDAnalysis/mdacli/actions/workflows/py38.ym

 .. |test_py39| image:: https://github.com/MDAnalysis/mdacli/actions/workflows/py39.yml/badge.svg
   :alt: Github Actions Build Status
   :target: https://github.com/MDAnalysis/mdacli/actions/workflows/py39.yml