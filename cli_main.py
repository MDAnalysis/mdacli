"""
Main entry point for the MDAnalysis CLI interface.

This also demonstrates how other third party libraries could incorporate
this functionality.
"""
# NOTE: names in this file are orientative
import argparse
import importlib
import inspect
import sys
import warnings

import MDAnalysis as mda
from MDAnalysis.analysis import __all__
from MDAnalysis.analysis.base import AnalysisBase
from numpydoc.docscrape import NumpyDocString


# modules in MDAnalysis.analysis packages that are ignored by MDA-CLI
# relevant modules used in this CLI factory
skip_mods = ('base', 'rdf_s')
relevant_modules = (_mod for _mod in __all__ if _mod not in skip_mods)

# global dictionary storing the parameters for all Analysis classes
analysis_interfaces = {}

# serves CLI factory
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


def add_to_CLIs(callable_obj, storage_dict):
    """Inspect Analysis class or function."""
    storage_dict[callable_obj.__name__] = {}
    storage_dict[callable_obj.__name__]["callable"] = callable_obj

    sig = inspect.signature(callable_obj)
    doc = NumpyDocString(callable_obj.__doc__)

    # args for CLIs
    positional_args = {}
    optional_args = {}

    for sig_name, sig_param in sig.parameters.items():

        if sig_param.kind == inspect.Parameter.VAR_KEYWORD:
            # this is **kwargs
            # define behaviour here
            # i doubt this is of any use,
            # if exists in analysis classes mostlikely refers to
            # receiving parameters remaining from other contexts
            pass

        elif sig_param.kind == inspect.Parameter.VAR_POSITIONAL:
            # this is for *args
            # I think the same rationale as for VAR_KEYWORD applies
            pass

        elif sig_param.default == inspect.Parameter.empty:
            for doc_param in doc["Parameters"]:
                if doc_param.name == sig_name:
                    positional_args[sig_name] = {
                        "type": doc_param.type.split()[0],
                        "desc": " ".join(doc_param.desc),
                    }
                    break
            # else reaches if the parameter in the signature is not present in the docstring
            # it shouldn't, but just in case :-)
            # unless we explicitly decide not to consider any parameters not referenced in the
            # documentation.
            else:
                # str is the default value of argparse arguments type parameter
                positional_args[sig_name] = {
                    "type": "str",
                    "desc": "No description available.",
                }

        else:
            for doc_param in doc["Parameters"]:
                if doc_param.name == sig_name:
                    optional_args[sig_name] = {
                        "type": doc_param.type.split()[0],
                        "default": sig_param.default,
                        "desc": " ".join(doc_param.desc),
                    }
                    break
            else:
                optional_args[sig_name] = {
                    "type": type(sig_param.default).__name__,  # corrected here
                    "default": sig_param.default,
                    "desc": "No description available.",
                }

    storage_dict[callable_obj.__name__]["positional"] = positional_args
    storage_dict[callable_obj.__name__]["optional"] = optional_args
    storage_dict[callable_obj.__name__]["desc"] = doc["Summary"]
    storage_dict[callable_obj.__name__]["desc_long"] = doc["Extended Summary"]

    # we can add here whatever we need more
    return callable_obj


def add_interface_CLI(cli_parser, interface_name, parameters):
    """
    Add subparsers to `cli_parser`.

    Parameters
    ----------
    cli_parser : argparse.sub_parser
        The main parser where the new parser will be added.

    interface_name : str
        Name of the interface name.

    parameters : dict
        Parameters needed to fill the argparse requirements for the
        CLI interface.
    """
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
    return


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


def maincli(ap):
    """Execute main client interface."""
    if len(sys.argv) < 2:
        ap.error("A subcommand is required.")

    try:
        args = ap.parse_args()
        main(**vars(args))
    except Exception as e:
        sys.exit("{}Error: {}{}".format(bcolors.fail, e, bcolors.endc))


ap = argparse.ArgumentParser()
cli_parser = ap.add_subparsers(title="MDAnalysis Analysis CLI")

# populates analysis_interfaces dictionary
for module in relevant_modules:
    module = importlib.import_module('MDAnalysis.analysis.' + module)
    for name, member in inspect.getmembers(module):
        if inspect.isclass(member) and issubclass(member, AnalysisBase):
            add_to_CLIs(member, analysis_interfaces)

# adds each Analysis class/function as a CLI under 'cli_parser'
# to be writen
for interface_name, parameters in analysis_interfaces.items():
    add_interface_CLI(cli_parser, interface_name, parameters)

# the entry point for this file needs to be added also to the
# setup.py file
if __name__ == "__main__":
    maincli(ap)
