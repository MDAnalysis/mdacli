#!/usr/bin/env python3
# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
#
# Copyright (c) 2021 Authors and contributors
#
# Released under the GNU Public Licence, v2 or any higher version
# SPDX-License-Identifier: GPL-2.0-or-later
"""Manage data saving."""
import json
import os
import sys
import zipfile

import numpy as np

from MDAnalysis.analysis.base import Results

from mdacli.colors import Emphasise


def save_results(fprefix, results):
    """
    Save the attributes of a results instance to disk.

    1D, 2D and 3D numpy arrays are saved to csv files. 1D arrays of the
    same length are vertically stacked to create a table. 2D arrays are
    saved directly. 3D arrays are split into 2D arrays along the
    shortest dimension and one CSV is saved for each 2D array created,
    and resulting CSVs are stored together in a ZIP file. Note: higher
    dimensional arrays are ignored.

    We try to save everything else in a JSON file. Non-serializable
    types are ignored.

    Parameters
    ----------
    fprefix : str
        prefix for all files saved

    results : `~MDAnalysis.analysis.base.Results`
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
            n_dims = item.ndim

            if n_dims == 1:
                list_1D.append(item)
                list_1D_labels.append(key)

            elif n_dims == 2:
                np.savetxt(
                    fname=f"{fprefix}_{key}.csv",
                    X=item,
                    fmt='%-10s',
                    #delimiter='\t',
                    )

            elif n_dims == 3:
                min_dim = np.argmin(item)
                files_to_zip = []

                # Split array along the dimension with smallest number
                # of entries
                splitted_item = np.split(
                    item,
                    item.shape[min_dim],
                    axis=min_dim,
                    )

                for arr in splitted_item:
                    files_to_zip.append(f"{key}_dim_{min_dim}_idx_{i}.csv")
                    np.savetxt(
                        fname=files_to_zip[-1],
                        X=np.squeeze(arr),
                        delimiter=',',
                        )

                # Compress all csv files into a single zip archive
                with zipfile.ZipFile(f'{key}.zip', 'w') as zipF:
                    for file_name in files_to_zip:
                        zipF.write(
                            file_name,
                            compress_type=zipfile.ZIP_DEFLATED,
                            )
                        os.remove(file_name)

            elif n_dims > 3:
                np.save(
                    f"{fprefix}_{key}",
                    item,
                    allow_pickle=False,
                    )

            else:
                warnings.warn(
                    Emphasise.warning(
                        "Saving numpy arrays with more than "
                        "three dimensions is currently not supported."
                        )
                    )

        elif (isinstance(item, (bool, int, float, list, tuple, dict))
              or item is None):
            # This can be encoded in a json file
            json_dict[key] = item

        else:
            warnings.warn(Emphasise.warning(f"Saving {key} of type {type(item)}"
                          "is currently not supported."))

    # Stack 1D arrays and save teheem to csv
    if len(list_1D) > 0:
        out_lists, out_lables = stack_1d_arrays_list(list_1D, list_1D_labels)

        command_run = ' '.join(sys.argv) + "\n"

        for out_list, out_label in zip(out_lists, out_lables):
            out_label = out_label.flatten()

            # [3:] to align lables with entries
            savetxt_w_command(
                fname=f"{fprefix}_{'_'.join(out_label)}.csv",
                X=out_list.T,
                header=''.join([f"{i:>25}" for i in out_label])[3:]
                )

    # Save everything which is left to a json file
    # if json dict is not empty
    if json_dict:
        with open(f'{fprefix}.json', 'w') as f:
            json.dump(json_dict, f)

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




def savetxt_w_command(fname, X, header='', fsuffix=".csv", **kwargs):
    """
    Save CSV data with info about execution command.

    Adds the command line input to the header and checks for a doubled
    defined filesuffix.
    """
    header = "{}\n{}".format(get_cli_input(), header)
    fname = "{}{}".format(fname, (not fname.endswith(fsuffix)) * fsuffix)
    np.savetxt(fname, X, header=header, **kwargs) 


def get_cli_input():
    """Return a proper fomatted string of the command line input."""
    program_name = os.path.basename(sys.argv[0])
    # Add additional quotes for connected arguments.
    arguments = [
        '"{}"'.format(arg).strip()
        if " " in arg else arg
        for arg in sys.argv[1:]
        ]

    return "Command line was: {} {}".format(program_name, " ".join(arguments))
