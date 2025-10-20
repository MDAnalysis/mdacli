=======================
Philosophy and approach
=======================

The mdacli project evolved from our experiences in
`taurenmd <https://github.com/joaomcteixeira/taurenmd>`_ and
`maicos <https://gitlab.com/maicos-devel/maicos>`_ which
both build a CLI interface on the fly.

They were developed since
we were facing the same problem in the labs regularly.
New students/scientists start writing their analysis
scripts and face the same challenges and problems i.e.

* How to initialize the universe and loop through frames without copying many lines of code?
* How to write a CLI parser to analyze several of their simulations?
* How to process and save their trajectories in a clever way?
* ...

Some of these problems can be solved by using the
:class:`MDAnalysis.analysis.base.AnalysisBase` class of MDA. However,
this class is limited to python, and sometimes a direct CLI to these
scripts is very helpful for the day-to-day analysis.
A generic CLI wrapper for all classes based on the AnalysisBase could
therefore help people to analyze their simulation
data with the least effort. With this approach, it is easier to use for
MDA-users since they just stay within their known universe
with known selection commands and results structures.
An existing framework makes it also more attractive for users
and developers to write their analysis using the
`base.AnalysisBase`.

Starting from `taurenmd` and `maicos` we
developed a general CLI for any
:class:`MDAnalysis.analysis.base.AnalysisBase` class.
mdacli detects all analysis classes located inside the
MDA project and builds a CLI wrapper around them. The wrapper
is generic so it also applies to any downstream
project that uses the :class:`MDAnalysis.analysis.base.AnalysisBase`
as parent class. If new classes are added to the
MDA codebase they are just they will show up without
any adjustments to `mdacli` itself.

The core of the wrapper is a docstring parser in combination
with an argument inspection using the `inspect` library. Based on
a created dictionary containing
each parameter of the class with its docstring and type, the actual command
line interface is build using :py:mod:`argparse`. The syntax of the
topology and the trajectory flags (`-f`, `-s`, ...) is inspired by the
`GROMACS CLI <https://manual.gromacs.org/documentation/current/user-guide/cmdline.html>`_
syntax.

The interface also provides a
way to save the data using that all analysis results of an AnalysisBase
class is stored inside results objects. The saving routines automatically
detects the type of the results and saves them either as JSON dumps (for
simple variables), CSV files (for 1D and 2D arrays), or zipped data dumps
for high higher-dimensional arrays.
