=====
Usage
=====

To use ``mdacli``, after installation
open your terminal and run::

 mda -h

This command will provide a list of all available modules. A list of the
:ref:`Available modules` is also available on the page of the documentation. Ask
help `-h` in each module available for detailed instruction on how to use each
module command-line client.

`mdacli` modules' parameters emulate the parameters of the `Analysis` classes
from `MDAnalysis`. So, each module will have its own requirements. Some will
require two trajectories, others `AtomGroup` selections, etc. You will see that
all is explained in each client `-h` option.

For example, to calculate a radial distribution function (RDF) between two groups use for
example::

 mda interrdf -h

A sample water trajectory is provided of rigid SPC/E water is
provided `online`_.
`topol.tpr` contains a GROMACS topolgy file and `traj.trr` is
the corresponding trajectory. The oxygen-oxygen
rdf can be calculated using::

 mda interrdf -s topol.tpr -f traj.trr -g1 "name OW" -g2 "name OW"

The oxygen atoms are selected with the `-g1` and `-g2` flags.

A more verbose output is achieved by using the `-v` flag. Even more
information is provided with the `--debug` flag.
All warnings
of the run can also directly stored to a log file using the `--logfile`
flag::

 mda --debug --logfile rdf.log interrdf -v -s topol.tpr -f traj.trr -g1 "name OW" -g2 "name OW"

The results of each analysis are stored by default in the current directory.
The output
directory can be changed with the `-o` flag and an additional prefix can be
set using the `-pre` flag.

The results of the RDF calculations
are two `.csv` files. The actual RDF is saved in the 2nd and 3rd columns
of `InterRDF_count_bins_rdf.csv`. The header rows of each file provide
information about the stored data. Simple results such as bare numbers or
strings are stored as `JSON` dumps. More complex data such as
4 or higher-dimensional arrays are saved as a bunch of CSV files zipped
together. A similar procedure will happen for each module.

.. _online: https://github.com/MDAnalysis/mdacli/tree/main/data

