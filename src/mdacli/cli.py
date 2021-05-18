#!/usr/bin/env python3
# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
#
# Copyright (c) 2021 Authors and contributors
#
# Released under the GNU Public Licence, v2 or any higher version
# SPDX-License-Identifier: GPL-2.0-or-later

"""
Main entry point for the MDAnalysis CLI interface.

This also demonstrates how other third party libraries could incorporate
this functionality.
"""
import argparse
import importlib
import inspect
import json
import os
import sys
import warnings
import zipfile

import MDAnalysis as mda
import numpy as np
from MDAnalysis.analysis import __all__
from MDAnalysis.analysis.base import AnalysisBase, Results

from mdacli.colors import Emphasise
from mdacli.utils import convert_str_time, parse_callable_signature, parse_docs


# modules in MDAnalysis.analysis packages that are ignored by mdacli
# relevant modules used in this CLI factory
# hydro* are removed here because they have a different folder/file structure
# and need to be investigated separately
skip_mods = ('base', 'hydrogenbonds', 'hbonds')
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


def _warning(message, *args, **kwargs):
    print(Emphasise.warning(f"Warning: {message}"))


warnings.showwarning = _warning


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

    # adds analyze_data function as the default func parameter.
    # this is possible because the analyze_data function is equal to all
    # Analysis Classes
    common_group.set_defaults(func=analyze_data)

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
        type=str,
        default="0",
        help="frame or start time for evaluation. (default: %(default)s)"
        )

    common_group.add_argument(
        "-e",
        dest="end",
        type=str,
        default="-1",
        help="frame or end time for evaluation. (default: %(default)s)"
        )

    common_group.add_argument(
        "-dt",
        dest="dt",
        type=str,
        default="1",
        help="step or time step for evaluation. (default: %(default)s)"
        )

    common_group.add_argument(
        "-pre",
        dest="output_prefix",
        type=str,
        default="",
        help="Additional prefix for all output files. Files will be "
             " automatically named by the used module (default: %(default)s)"
        )

    common_group.add_argument(
        "-o",
        dest="output_directory",
        type=str,
        default=".",
        help="Directory in which the output files produced will be stored."
             "(default: %(default)s)"
        )

    common_group.add_argument(
        "-v",
        dest="verbose",
        help="Be loud and noisy",
        action="store_true"
        )

    pos_ = sorted(list(parameters["positional"].items()), key=lambda x: x[0])
    opt_ = sorted(list(parameters["optional"].items()), key=lambda x: x[0])

    parameters_to_parse = pos_ + opt_
    mandatory_parameters_group = analysis_class_parser.add_argument_group(
        title="Mandatory Parameters",
        description="Mandatory parameters of this Analysis",
        )
    groups = len(pos_) * [mandatory_parameters_group]

    # Only create parser if optional arguments exist
    if len(opt_) > 0:
        optional_parameters_group = analysis_class_parser.add_argument_group(
            title="Optional Parameters",
            description="Optional parameters specific of this Analysis",
            )
        groups += len(opt_) * [optional_parameters_group]

    action_dict = {True: "store_false", False: "store_true"}
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


def stack_1d_arrays_list(list_1D, extra_list=None):
    """Stack a list of 1D numpy arrays of the same length vertically together.

    The result is a list containing 2D arrays where each array got the same
    number of rows.

    Parameters
    ----------
    list_1d : list
        list of 1 dimensional numpy arrays

    extra_list : list
        additional list of numpy arrays on which the
        operations are executed as for ``list_1d``

    Returns
    -------
    out_list : list
        list of stacked 2D numpy arrays organized by their length
    out_extra : list
        list of stacked 2D numpy applied applied to the same operations
        as out_list
    """
    # Sort for lengths
    lengths = np.array([len(a) for a in list_1D])
    sorted_idx = np.argsort(lengths)

    # Sort lists according to the lengths of the items
    list_1D_sorted = [list_1D[i] for i in sorted_idx]

    # Count the number of items for each length
    counts = np.unique(lengths, return_counts=True)[1]
    new_length_idx = np.hstack([[0], np.cumsum(counts)])

    out_lists = []
    # Concentanate lists of the same lenngth
    for i in range(0, len(new_length_idx) - 1):
        out_lists.append(np.vstack(list_1D_sorted[new_length_idx[i]:
                                                  new_length_idx[i + 1]]))

    if extra_list is not None:
        extra_list_sorted = [extra_list[i] for i in sorted_idx]
        out_extra = []
        for i in range(0, len(new_length_idx) - 1):
            out_extra.append(np.vstack(extra_list_sorted[new_length_idx[i]:
                                                         new_length_idx[i + 1]]
                                       ))

        return out_lists, out_extra
    else:
        return out_lists


def save_results(fprefix, results):
    """Save the attributes of a results instance to disk.

    1D, 2D and 3D numpy arrays are saved to csv files. All 1D arrays
    of the same lengths are veertically stacked. For 3D arrays
    a csv file is created for the dimension with the lowest number of
    indices. Higher dimensional arrays are ignored.

    Everything else is tried to saved inside a json file. Types which
    can not be saved into json are ignored.

    Parameters
    ----------
    fprefix : str
        prefix for all files saved
    results : `MDAnalysis.analysis.base.Results`
        A Results instance from which the stored data is taken.
    """
    list_1D = []
    list_1D_labels = []
    json_dict = {}

    for key, item in results.items():
        if isinstance(item, Results):
            # Run `save_results` recursively if
            # `item` is results instancee
            save_results(f"{fprefix}_{key}", item)
        elif isinstance(item, np.ndarray):
            # Remove extra dimensions
            item = np.squeeze(item)
            n_dims = len(item.shape)

            if n_dims == 1:
                list_1D.append(item)
                list_1D_labels.append(key)
            elif n_dims == 2:
                np.savetxt(fname=f"{fprefix}_{key}.csv",
                           X=item,
                           delimiter=',')
            elif n_dims == 3:
                min_dim = np.argmin(item)
                files_to_zip = []
                # Split array along the dimension with smallest number
                # of entries
                for i, arr in enumerate(np.split(item,
                                                 item.shape[min_dim],
                                                 axis=min_dim)):
                    files_to_zip.append(f"{key}_dim_{min_dim}_idx_{i}.csv")
                    np.savetxt(fname=files_to_zip[i],
                               X=np.squeeze(arr),
                               delimiter=',')

                # Compress all csv files into a single zip archive
                with zipfile.ZipFile(f'{key}.zip', 'w') as zipF:
                    for file in files_to_zip:
                        zipF.write(file, compress_type=zipfile.ZIP_DEFLATED)
                        os.remove(file)

            else:
                warnings.warn("Saving numpy arrays with more than "
                              "three dimensions is currently not supported.")
        elif (isinstance(item, (bool, int, float, list, tuple, dict))
              or item is None):
            # This can be encoded in a json file
            json_dict[key] = item

        else:
            warnings.warn(f"Saving {key} of type {type(item)}"
                          "is currently not supported.")

    # Stack 1D arrays and save teheem to csv
    if len(list_1D) > 0:
        out_lists, out_lables = stack_1d_arrays_list(list_1D, list_1D_labels)

        for out_list, out_label in zip(out_lists, out_lables):
            out_label = np.squeeze(out_label).tolist()

            # [3:] to align lables with entries
            np.savetxt(fname=f"{fprefix}_{'_'.join(out_label)}.csv",
                       X=out_list.T,
                       header=''.join([f"{i:>25}" for i in out_label])[3:]
                       )

    # Save everything which is left to a json file
    with open(f'{fprefix}.json', 'w') as f:
        json.dump(json_dict, f)


def analyze_data(
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
    """Perform main client logic."""
    u = mda.Universe(topology, trajectories)

    # so, here we need to do some investigation, questions are: * do all
    # Analysis classes/functions have the same execution interface?  * we need
    # to discuss with Oliver maybe to ensure that the above point is true *
    # otherwise we would need to write a dedicated cli main function to each
    # Analysis class * however, we can accept a certain number of different
    # interfaces and we can handle the control flow with try/catch statements
    # until the correct interface is found. try/catch on polymorphism, yeah! :-D

    # Convert special types (i.e AtomGroups)
    # Ugly that we have to parse again... but currently I have no better idea :(
    params = parse_docs(analysis_callable)[2]  # Index [2] for paramaters
    for param_name, dictionary in params.items():
        if "AtomGroup" in dictionary['type']:
            sel = u.select_atoms(analysis_kwargs[param_name])
            if len(sel) > 0:
                analysis_kwargs[param_name] = sel
            else:
                raise ValueError(
                    "AtomGroup `-{}` with selection `{}` does not "
                    "contain any atoms".format(
                        param_name,
                        analysis_kwargs[param_name]
                        )
                    )
        elif "Universe" in dictionary['type']:
            analysis_kwargs[param_name] = u

    with warnings.catch_warnings():
        warnings.simplefilter('always')
        startframe = convert_str_time(analysis_kwargs.pop("begin"), u.trajectory.dt)  # noqa: E501
        stopframe = convert_str_time(analysis_kwargs.pop("end"), u.trajectory.dt)  # noqa: E501
        step = convert_str_time(analysis_kwargs.pop("dt"), u.trajectory.dt)

        # raises error if frame selection range is an empty selection
        if not list(range(u.trajectory.n_frames)[slice(startframe, stopframe, step)]):  # noqa: E501
            raise ValueError("Trajectory frame range {}:{}:{} is not valid for {} frames."  # noqa: E501
                             "".format(startframe, stopframe, step, u.trajectory.n_frames))  # noqa: E501

    # Collect paramaters not necessary for initilizing ac object.
    analysis_kwargs.pop("func")
    verbose = analysis_kwargs.pop("verbose")
    output_directory = analysis_kwargs.pop("output_directory")
    output_prefix = analysis_kwargs.pop("output_prefix")
    output_prefix += "_" if len(output_prefix) > 0 else ""

    ac = analysis_callable(**analysis_kwargs)
    ac.run(start=startframe,
           stop=stopframe,
           step=step,
           verbose=verbose)

    save_results(os.path.join(output_directory,
                              f"{output_prefix}{type(ac).__name__}"),
                 ac.results)


def maincli(ap):
    """Execute main client interface."""
    if len(sys.argv) < 2:
        ap.error("A subcommand is required.")

    try:
        args = ap.parse_args()
        analyze_data(**vars(args))
    except Exception as e:
        sys.exit(Emphasise.error(f"Error: {e}"))


def setup_clients():
    """
    Set up ArgumentParser clients.

    Returns
    -------
    argparse.ArgumentParser instance
    """
    ap = argparse.ArgumentParser()
    cli_parser = ap.add_subparsers(title="MDAnalysis Analysis CLI")

    # populates analysis_interfaces dictionary
    for module in relevant_modules:
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            module = importlib.import_module('MDAnalysis.analysis.' + module)
        for _, member in inspect.getmembers(module):
            if inspect.isclass(member) and issubclass(member, AnalysisBase) \
               and member is not AnalysisBase:
                parse_callable_signature(member, analysis_interfaces)

    # adds each Analysis class/function as a CLI under 'cli_parser'
    # to be writen
    for interface_name, parameters in analysis_interfaces.items():
        create_CLI(cli_parser, interface_name, parameters)

    return ap


def main():
    """Execute main CLI entry point."""
    maincli(setup_clients())
