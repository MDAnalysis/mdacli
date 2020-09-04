"""
Main entry point for the MDAnalysis CLI interface.

This also demonstrates how other third party libraries could incorporate
this functionality.
"""
import argparse
import importlib
import inspect
import re
import sys
import warnings
from collections import defaultdict

import MDAnalysis as mda
from MDAnalysis.analysis import __all__
from MDAnalysis.analysis.base import AnalysisBase


# modules in MDAnalysis.analysis packages that are ignored by MDA-CLI
# relevant modules used in this CLI factory
# hydro* are removed here because they have a different folder/file structure
# and need to be investigated separately
skip_mods = ('base', 'rdf_s', 'hydrogenbonds', 'hbonds')
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


def parse_callable_signature(callable_obj, storage_dict):
    """
    Parse a callable object to a convenient dictionary for CLI creation.

    The parameters used in the CLI are a combination of the callable
    signature and the information in the callable docstring.

    Parameters
    ----------
    callable_obj : callable
        The callable object to inspect. Details of this object required
        for the creation of a CLI are added to the `storage_dict`.

    storage_dict : dict
        The dictionary that stores the details of the callable.

    Returns
    -------
    None
        Modifies `storage_dict` in place.
    """
    storage_dict[callable_obj.__name__] = {}
    storage_dict[callable_obj.__name__]["callable"] = callable_obj

    sig = inspect.signature(callable_obj)
    summary, summary_extended, doc = parse_docs(callable_obj)

    # args for CLIs
    positional_args = {}
    optional_args = {}

    # tuple to ignore args and kwargs, these are note necessary to consider
    # for the purpose of this implementation - at least to currently
    ARGS_KWARGS = (
        inspect.Parameter.VAR_KEYWORD,
        inspect.Parameter.VAR_POSITIONAL,
        )

    # brings to local scope
    EMPTY = inspect.Parameter.empty

    # for each parameter in the callable's signature
    for sig_name, sig_param in sig.parameters.items():

        if sig_param.kind in ARGS_KWARGS:
            pass

        # positional parameter
        elif sig_param.default == EMPTY:

            # parameter type and description is extract from docstring
            for param_name, doc_param in doc.items():
                if param_name == sig_name:
                    positional_args[sig_name] = {
                        "type": doc_param['type'],
                        "desc": doc_param['desc'],
                        }
                    break  # done, jumps off the loop

            else:
                # else reaches if the parameter in the signature is not present in
                # the docstring. It shouldn't, but just in case :-)
                # unless we explicitly decide not to consider any parameters not
                # referenced in the documentation, this should be kept
                #
                # str is the default value of argparse arguments type parameter
                positional_args[sig_name] = {
                    "type": "str",
                    "desc": "No description available.",
                    }

        # named parameters
        else:
            for param_name, doc_param in doc.items():
                if param_name == sig_name:
                    optional_args[sig_name] = {
                        "type": doc_param['type'],  # type taken form docstring
                        "default": sig_param.default,  # but default taken from signature
                        "desc": doc_param['desc'],
                    }
                    break  # parameter captured, break the loop
            else:
                # if the parameter is in signature but NOT in the docstring
                # uses type.__name__ to match with STR_TYPE_DICT keys
                optional_args[sig_name] = {
                    "type": type(sig_param.default).__name__,  # corrected here
                    "default": sig_param.default,
                    "desc": "No description available.",
                    }

    # places all information captured for the callable in the dictionary
    storage_dict[callable_obj.__name__]["positional"] = positional_args
    storage_dict[callable_obj.__name__]["optional"] = optional_args
    storage_dict[callable_obj.__name__]["desc"] = summary
    storage_dict[callable_obj.__name__]["desc_long"] = summary_extended

    return


def parse_docs(klass):
    """
    Parse classes docstrings to a convenient dictionary.

    This parser is based on NumpyDocString format, yet it is not so
    strict. Combined docstrings from class main docstring and `__init__`
    method.

    Parameters
    ----------
    klass : callable
        A klass object from which a DOCSTRING can be extracted.

    Returns
    -------
    tuple (str, str, dict of dict)
        * One line summary description of the callable
        * Extended description of the callable
        * dictionary where keys are parameter names and subdictinary
            has keys "type" and "desc" for parameter type and description.
    """
    doc = klass.__doc__ or ''
    doc += klass.__init__.__doc__ or ''

    doc_lines = [s for s in (s.strip() for s in doc.lstrip().split('\n')) if s]

    # first docstring sentence is the summary
    summary = doc_lines[0]

    # sometimes signature parameters in docstring are referred as "Arguments"
    try:
        param_index = doc_lines.index('Parameters')
    except ValueError:
        param_index = doc_lines.index('Arguments')

    # the extended summary is every text that exists between the summary
    # and the Parameters title line
    summary_extended = '\n'.join(doc_lines[1:param_index])

    # the line to start collecting parameters is `param_index` + 2
    # because of the "Parameters" title and the '---------' underscore
    par_i = param_index + 2

    # search for the line where Parameters section ends
    # will search until the end of the docstring or until a '----'-like
    # line is found - corresponding to the start of another section
    for i, line in enumerate(doc_lines[par_i:], start=par_i):
        if '----' in line:  # at least of Note\n----
            # -1 because the exact line is the one of the title not of
            # the underlines
            end_param_line = i - 1
            break
    else:
        # if the end of the docstring is reached, takes the last line
        end_param_line = -1

    # starts collecting parameters in doctrings
    params = defaultdict(dict)  # the dictionary which will be returned

    # temporary parameter description list
    desc_tmp = []

    # regex to find parameter types
    type_regex = re.compile(r'^(\w+|\{.*\})')

    # goes back to front to register descriptions first ;-)
    # considers only the Parameters section
    for line in doc_lines[par_i: end_param_line][::-1]:
        if ' : ' in line:
            par_name, others_ = line.split(' : ')
            par_type = type_regex.findall(others_)[0]
            params[par_name]['type'] = par_type
            params[par_name]['desc'] = ' '.join(desc_tmp[::-1])
            desc_tmp.clear()
        else:
            desc_tmp.append(line)

    return summary, summary_extended, params


def create_CLI(cli_parser, interface_name, parameters):
    """
    Add subparsers to `cli_parser`.

    Subparsers parameters are divided in three categories:

    1) Common Analysis classes Parameters
        Common to all generated CLIs, for example:
            * topology, trajectory, time frame

    2) Mandatory Parameters
        mandatory parameters are defined in the CLI as named parameters
        as per design

    3) Optional Parameters
        Named parameters in the Analysis class

    All CLI's parameters are named parameters.

    Parameters
    ----------
    cli_parser : argparse.sub_parser
        The main parser where the new parser will be added.

    interface_name : str
        Name of the interface name.

    parameters : dict
        Parameters needed to fill the argparse requirements for the
        CLI interface.

    Returns
    -------
    None
    """
    # creates the subparser
    analysis_class_parser = cli_parser.add_parser(
        interface_name,
        help=parameters["desc"],
        description=parameters["desc"] + "\n\n" + parameters["desc_long"],
        formatter_class=argparse.RawDescriptionHelpFormatter,
        )

    common_group = analysis_class_parser.add_argument_group(
        title="Common Analysis Parameters",
        )

    # adds main function as the default func parameter.
    # this is possible because the main function is equal to all Analysis Classes
    common_group.set_defaults(func=main)

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
        # these can be added on the create_CLI level
        topology,
        trajectories,
        # analysis_callable paramter needs to be injected where from the
        # global dictionary using argparse.set_defaults()
        # https://docs.python.org/3/library/argparse.html#argparse.ArgumentParser.set_defaults
        analysis_callable=None,
        **analysis_kwargs,
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
    params = parse_docs(analysis_callable)[2]  # Index [2] for paramaters
    for param_name, dictionary in params.items():
        if "AtomGroup" in dictionary['type']:
            analysis_kwargs[param_name] = u.select_atoms(analysis_kwargs[param_name])

    with warnings.catch_warnings():
        warnings.simplefilter('always')
        begin = analysis_kwargs.pop("begin")
        if begin > u.trajectory.totaltime:
            raise ValueError("Start ({:.2f} ps) is larer than total time "
                             "({:.2f} ps).".format(analysis_kwargs["begin"],
                                                   u.trajectory.totaltime))
        elif begin > 0:
            startframe = int(analysis_kwargs["begin"] // u.trajectory.dt)
        else:
            startframe = 0
    
        end = analysis_kwargs.pop("end")
        if end is not None:
            stopframe = int(analysis_kwargs["end"] // u.trajectory.dt)
            analysis_kwargs["end"] += 1  # catch also last frame in loops
        else:
            stopframe = None

        dt = analysis_kwargs.pop("dt")
        if dt > 0:
            step = int(analysis_kwargs["dt"] // u.trajectory.dt)
        else:
            step = 1

    # Collect paramaters not necessary for initilizing ac object.
    verbose = analysis_kwargs.pop("verbose")
    analysis_kwargs.pop("func")

    ac = analysis_callable(**analysis_kwargs)
    results = ac.run(start=startframe,
                     stop=stopframe,
                     step=step,
                     verbose=verbose)

    # prototype lines to test functionality TO REMOVE
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


def setup_clients():
    """
    Setup ArgumentParser clients.

    Returns
    -------
    argparse.ArgumentParser instance
    """
    ap = argparse.ArgumentParser()
    cli_parser = ap.add_subparsers(title="MDAnalysis Analysis CLI")

    # populates analysis_interfaces dictionary
    for module in relevant_modules:
        module = importlib.import_module('MDAnalysis.analysis.' + module)
        for name, member in inspect.getmembers(module):
            if inspect.isclass(member) and issubclass(member, AnalysisBase):
                parse_callable_signature(member, analysis_interfaces)

    # adds each Analysis class/function as a CLI under 'cli_parser'
    # to be writen
    for interface_name, parameters in analysis_interfaces.items():
        create_CLI(cli_parser, interface_name, parameters)

    return ap

# the entry point for this file needs to be added also to the
# setup.py file
if __name__ == "__main__":
    maincli(setup_clients())
