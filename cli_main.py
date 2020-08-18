"""
Main entry point for the MDAnalysis CLI interface.

This also demonstrates how other third party libraries could incorporate
this functionality.
"""

# NOTE: names in this file are orientative
import argparse
import sys
import warnings
import importlib
import inspect
from pprint import pprint

import MDAnalysis as mda
from MDAnalysis.analysis import __all__
from MDAnalysis.analysis import analysis_interfaces
from MDAnalysis.analysis.base import add_to_CLIs, AnalysisBase


skip_mods = ('base', 'rdf_s')
relevant_modules = (_mod for _mod in __all__ if _mod not in skip_mods)

for module in relevant_modules:
    module = importlib.import_module('MDAnalysis.analysis.' + module)
    for name, member in inspect.getmembers(module):
        if inspect.isclass(member) and issubclass(member, AnalysisBase):
            add_to_CLIs(member)


# Coloring for warnings and errors
class bcolors:
    warning = '\033[93m'
    fail = '\033[91m'
    endc = '\033[0m'


def _warning(message,
             category=UserWarning,
             filename='',
             lineno=-1,
             file=None,
             line=None):
    print("{}Warning: {}{}".format(bcolors.warning, message, bcolors.endc))


warnings.showwarning = _warning


#@add_to_CLIs
#class NewAnalysis(AnalysisBase):
#    """
#    This is a doc.
#
#    The trajectory is read, frame by frame, and the atoms in `atomgroup` are
#    histogrammed on a 3D grid with spacing `delta`.
#
#    Parameters
#    ----------
#    atomgroup : AtomGroup or UpdatingAtomGroup
#            Group of atoms (such as all the water oxygen atoms) being analyzed.
#            This can be an :class:`~MDAnalysis.core.groups.UpdatingAtomGroup` for
#            selections that change every time step.
#    atomgroup2 : AtomGroup
#            fancy second atomgroup
#    delta : float (optional)
#            Bin size for the density grid in ångström (same in x,y,z).
#    float_arg : list or str or tuple
#        Does something.
#    bool_arg : bool
#        a bool argument
#    """
#
#    def __init__(
#        self,
#        atomgroup,
#        atomgroup2,
#        par1,
#        par2,
#        par3,
#        par4,
#        *args,
#        float_arg=0.5,
#        none_arg=None,
#        bool_arg=True,
#        tuple_arg=(4, 5),
#        str_arg="string",
#        **kwargs,
#    ):
#        super().__init__(atomgroup.universe.trajectory)
#        self._par1 = par1
#        self._par2 = par2
#        self._par3 = par3
#        self._par4 = par4
#
#    def _single_frame(self):
#        pass
#
#
## testing another analysis
#@add_to_CLIs
#class NewAnalysis2(NewAnalysis):
#    """
#    This is another doc.
#
#    The trajectory is read, frame by frame, and the atoms in `atomgroup` are
#    histogrammed on a 3D grid with spacing `delta`.
#
#    Parameters
#    ----------
#    atomgroup : AtomGroup or UpdatingAtomGroup
#            Group of atoms (such as all the water oxygen atoms) being analyzed.
#            This can be an :class:`~MDAnalysis.core.groups.UpdatingAtomGroup` for
#            selections that change every time step.
#    atomgroup2 : AtomGroup
#            fancy second atomgroup
#    delta : float (optional)
#            Bin size for the density grid in ångström (same in x,y,z).
#    float_arg : list or str or tuple
#        Does something.
#    """
#
#    def _single_frame(self):
#        pass


###############################################################################

STR_TYPE_DICT = {
    "bool": bool,
    "str": str,
    "list": list,
    "tuple": tuple,
    "int": int,
    "float": float,
    "complex": complex,
    "NoneType": type(None),
    "AtomGroup": mda.AtomGroup,
}


def add_interface_CLI(cli_parser, interface_name, parameters):
    """Function for adding subparsers to cli_parser"""
    #pprint(parameters)

    analysis_class_parser = cli_parser.add_parser(
        interface_name, help="".join(parameters["desc"])
    )

    analysis_class_parser.description = " ".join(
        parameters["desc"] + parameters["desc_long"]
    )

    common_group = analysis_class_parser.add_argument_group(
        title="Common Analysis Parameters"
    )

    # adds main function as the default func parameter.
    # this is possible because the main function is equal to all Analysis Classes
    common_group.set_defaults(func=main)

    # Adds also the callable
    common_group.set_defaults(analysis_callable=parameters["callable"])

    common_group.add_argument(
        "-s",
        dest="topology",
        type=str,
        default="topol.tpr",
        help="The topolgy file. "
        "The FORMATs {} are implemented in MDAnalysis."
        "".format(", ".join(mda._PARSERS.keys())),
    )

    common_group.add_argument(
        "-f",
        dest="trajectories",
        type=str,
        default=None,
        nargs="+",
        help="A single or multiple trajectory files. "
        "The FORMATs {} are implemented in MDAnalysis."
        "".format(", ".join(mda._READERS.keys())),
    )

    common_group.add_argument(
        "-b",
        dest="begin",
        type=float,
        default=0,
        help="start time (ps) for evaluation. (default: %(default)s)"
    )

    common_group.add_argument(
        "-e", dest="end", type=float, default=None,
        help="end time (ps) for evaluation. (default: %(default)s)"
    )

    common_group.add_argument(
        "-dt",
        dest="dt",
        type=float,
        default=0,
        help="time step (ps) to read analysis frame. "
             "If `0` take all frames (default: %(default)s)"
    )

    common_group.add_argument(
        "-v",
        dest="verbose",
        help="Be loud and noisy",
        action="store_true")

    pos_ = sorted(list(parameters["positional"].items()), key=lambda x: x[0])
    opt_ = sorted(list(parameters["optional"].items()), key=lambda x: x[0])

    parameters_to_parse = pos_ + opt_

    mandatory_parameters_group = analysis_class_parser.add_argument_group(
        title="Mandatory Parameters",
        description="Mandatory parameters of this Analysis",
    )

    optional_parameters_group = analysis_class_parser.add_argument_group(
        title="Optional Parameters",
        description="Optional parameters specific of this Analysis",
    )

    groups = len(pos_) * [mandatory_parameters_group] + \
             len(opt_) * [optional_parameters_group]

    action_dict = {True: "store_false", False: "store_true"}
    #print(parameters_to_parse)
    for group, (name, args_dict) in zip(groups, parameters_to_parse):

        # prepares parameters before add_argument
        try:
            type_ = STR_TYPE_DICT[args_dict["type"]]
        except KeyError:
            type_ = str

        try:
            default = args_dict["default"]
        except KeyError:
            default = None

        name_par = "-" + name

        description = args_dict["desc"]

        if issubclass(type_, (list, tuple)):
            group.add_argument(
                name_par, dest=name, nargs="+", default=default,
                help="{} (default: %(default)s)".format(description)
            )
        elif type_ is bool:
            group.add_argument(
                name_par,
                dest=name,
                action=action_dict[default],
                default=default,
                help=description,
            )
        elif type_ is mda.AtomGroup:
            group.add_argument(
                name_par,
                dest=name,
                type=str,
                default=default,
                help=description + " Use a MDAnalysis selection string."
            )
        else:
            group.add_argument(
                name_par, dest=name, type=type_, default=default,
                help="{} (default: %(default)s)".format(description)
            )


def main(
    # top and trajs need to be positional parameters in all CLIs
    # these can be added on the add_interface_CLI level
    topology,
    trajectories,
    # analysis_callable paramter needs to be injected where from the
    # global dictionary using argparse.set_defaults()
    # https://docs.python.org/3/library/argparse.html#argparse.ArgumentParser.set_defaults
    analysis_callable=None,
    **analysis_kwargs
):
    """
    Main client logic.
    """
    u = mda.Universe(topology, trajectories)

    # so, here we need to do some investigation, questions are:
    # * do all Analysis classes/functions have the same execution interface?
    # * we need to discuss with Oliver maybe to ensure that the above point is true
    # * otherwise we would need to write a dedicated cli main function to each Analysis class
    # * however, we can accept a certain number of different interfaces and we can
    # handle the control flow with try/catch statements until the correct interface
    # is found. try/catch on polymorphism, yeah! :-D

    # Convert special types (i.e AtomGroups)
    # Ugly that we have to parse again... but currently I have no better idea :(
    for doc_param in NumpyDocString(analysis_callable.__doc__)["Parameters"]:
        if "AtomGroup" in doc_param.type:
            analysis_kwargs[doc_param.name] = u.select_atoms(analysis_kwargs[doc_param.name])

    ac = analysis_callable(**analysis_kwargs)

    with warnings.catch_warnings():
        warnings.simplefilter('always')
        if analysis_kwargs["begin"] > u.trajectory.totaltime:
            raise ValueError("Start ({:.2f} ps) is larer than total time "
                             "({:.2f} ps).".format(analysis_kwargs["begin"],
                                                   u.trajectory.totaltime))
        elif analysis_kwargs["begin"] > 0:
            startframe = int(analysis_kwargs["begin"] // u.trajectory.dt)
        else:
            startframe = 0
        if analysis_kwargs["end"] is not None:
            stopframe = int(analysis_kwargs["end"] // u.trajectory.dt)
            analysis_kwargs["end"] += 1  # catch also last frame in loops
        else:
            stopframe = None
        if analysis_kwargs["dt"] > 0:
            step = int(analysis_kwargs["dt"] // u.trajectory.dt)
        else:
            step = 1

    results = ac.run(start=startframe,
                     stop=stopframe,
                     step=step,
                     verbose=analysis_kwargs["verbose"])

    print(analysis_kwargs)
    sys.exit("Analysis complete. exiting...")
    # extract results?
    # here the same, how are the results collected?
    # for RMSD we need to access the 'rmsd' attribute after .run() method.
    # do we need to alter (add) a common interface on all the MDAnalysis
    # interfaces. This would imply a intervention on the MDA code itself.
    # we can definitively create a common method on all classes that links
    # to the classes specific method/attribute where the results are stored.

    save_results_to_some_file(results)


ap = argparse.ArgumentParser()
cli_parser = ap.add_subparsers(title="MDAnalysis Analysis CLI")


# adds each Analysis class/function as a CLI under 'cli_parser'
# to be writen
for interface_name, parameters in analysis_interfaces.items():
    add_interface_CLI(cli_parser, interface_name, parameters)


def maincli():
    """Execute main client interface."""
    if len(sys.argv) < 2:
        ap.error("A subcommand is required.")

    try:
        args = ap.parse_args()
        print(args)
        main(**vars(args))
    except Exception as e:
        sys.exit("{}Error: {}{}".format(bcolors.fail, e, bcolors.endc))


# the entry point for this file needs to be added also to the
# setup.py file
if __name__ == "__main__":
    maincli()
