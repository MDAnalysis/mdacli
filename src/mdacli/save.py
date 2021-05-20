#!/usr/bin/env python3
# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-

# Copyright (c) 2021 Authors and contributors
#
# Released under the GNU Public Licence, v2 or any higher version
# SPDX-License-Identifier: GPL-2.0-or-later
"""Manage data saving."""
import json
import os
import sys
import zipfile
from functools import partial
from pathlib import Path

import numpy as np
from MDAnalysis.analysis.base import Results


def save_results(results, fprefix="mdacli_results"):
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
    save_1D_arrays(results, fprefix=fprefix, remove=True)
    save_2D_arrays(results, fprefix=fprefix, remove=True)
    save_3D_arrays(results, fprefix=fprefix, remove=True)
    save_higher_dim_arrays(results, fprefix=fprefix, remove=True, min_ndim=4)
    save_json_serializables(results, remove=True, fname=fprefix)
    save_Results_object(results, fprefix=fprefix, remove=True)
    return


def save_1D_arrays(results, fprefix="1darray", remove=True):
    """
    Save 1D arrays from results.

    Parameters
    ----------
    results : dict-like
        Dictionary containing results.

    remove : bool
        If true remove keys mapping to 1D numpy arrays.
    """
    list_1D, list_1D_labels = get_1D_arrays(results)

    if not list_1D:
        return

    out_lists, out_lables = stack_1d_arrays_list(list_1D, list_1D_labels)

    for out_list, out_label in zip(out_lists, out_lables):
        out_label = out_label.flatten()

        # [3:] to align lables with entries
        savetxt_w_command(
            fname=f"{fprefix}_{'_'.join(out_label)}.csv",
            X=out_list.T,
            header=''.join([f"{i:>25}" for i in out_label])[3:]
            )

    return return_with_remove(results, list_1D_labels, remove)


def get_1D_arrays(results):
    """Get items from dict which correspond to np.ndarrays one dim."""
    list_1D = []
    list_1D_labels = []

    for key, value in results.items():
        if is_1d_array(value):
            list_1D.append(value)
            list_1D_labels.append(key)

    return list_1D, list_1D_labels


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


def save_2D_arrays(results, fprefix="2Darr", remove=True):
    """Save items of 2D array."""
    keys = []
    for key, value in results.items():
        value = try_to_squeeze_me(value)
        if is_2d_array(value):
            savetxt_w_command(fname=f"{fprefix}_{key}.csv", X=value)
            keys.append(key)

    return return_with_remove(results, keys, remove)


def save_3D_arrays(results, fprefix="3Darr", remove=True):
    """Save items of 2D array."""
    keys = []
    for key, value in results.items():
        value = try_to_squeeze_me(value)
        if is_3d_array(value):
            save_3D_array_to_2D_csv(
                value,
                arr_name=f"{fprefix}_{key}",
                zipit=True,
                )
            keys.append(key)

    return return_with_remove(results, keys, remove)


def save_3D_array_to_2D_csv(
        item,
        arr_name='arr',
        zipit=True,
        ):
    """
    Save 3D array to 2D CSVs.

    Has option to store all in a ZIP file.
    """
    min_dim = np.argmin(item)
    files_to_zip = []

    # Split array along the dimension with smallest number
    # of entries
    splitted_item = np.split(
        item,
        item.shape[min_dim],
        axis=min_dim,
        )

    save_to_folder = not zipit
    folder = None
    if save_to_folder:
        folder = Path(arr_name)
        folder.mkdir(parents=True, exist_ok=True)

    for i, arr in enumerate(splitted_item):
        fname = f"{arr_name}_dim_{min_dim}_idx_{i}.csv"
        foutname = folder.joinpath(fname) if folder else fname
        savetxt_w_command(fname=foutname, X=np.squeeze(arr))
        files_to_zip.append(fname)

    if zipit:
        save_files_to_zip(files_to_zip, zipname=arr_name, remove=True)

    return


def save_higher_dim_arrays(results, fprefix="XDarr", remove=True, min_ndim=4):
    """Save items of multidimensional arrays to CSV."""
    keys = []
    for key, value in results.items():
        value = try_to_squeeze_me(value)
        if is_higher_dimension_array(value, min_ndim):
            save_result_array(value, fprefix=fprefix, arr_name=key)
            keys.append(key)

    return return_with_remove(results, keys, remove)


def save_result_array(arr, fprefix='prefix'):
    """Save array to disk accoring to num of dimensions."""
    item = np.squeeze(arr)

    save_options = {
        item.ndim == 1: save_1D_arrays,
        item.ndim == 2: save_2D_arrays,
        item.ndim == 3: save_3D_arrays,
        item.ndim > 3: save_higher_dim_arrays,
        }

    save_options[True](item, fprefix=fprefix)
    return


def save_json_serializables(results, remove=True, **jsonargs):
    """Save serializable items to a JSON."""
    json_dict = {
        key: value
        for key, value in results.items()
        if is_serializable(value)
        }

    if remove:
        for key in json_dict.keys():
            results.pop(key)

    if json_dict:
        json_dict["command"] = get_cli_input()
        save_to_json(json_dict, **jsonargs)


def is_serializable(value):
    """Assert if value is json serializable."""
    try:
        json.dumps(value)
        return True
    except (TypeError, OverflowError):
        return False


def save_to_json(json_dict, fname='jdict', indent=4, sort_keys=True):
    """Save dictionary to JSON file."""
    with open(f'{fname}.json', 'w') as f:
        json.dump(json_dict, f, indent=indent, sort_keys=sort_keys)


def save_Results_object(results, fprefix='results', remove=True):
    """Save results if they are Results objects."""
    keys = []
    for key, value in results.items():
        if isinstance(value, Results):
            save_results(f"{fprefix}_{key}", value)
            keys.append(key)

    return return_with_remove(results, keys, remove)


def return_with_remove(ddict, keys, remove):
    """
    Serve all saving functions.

    If remove is true,
    Returns subset of keys from dict.
    Removes keys subset from original dict.

    Else, return None.
    """
    if remove:
        return {key: ddict.pop(key) for key in keys}

    else:
        return None


def save_files_to_zip(files, zipname='thezip', remove=True):
    """
    Compress all files into a single zip archive.

    Parameters
    ----------
    files : list-like
        File names to save to the ZIP archive.

    zipname : str
        The name of the zip file without extension.

    remove : bool, option, default True
        Removes the original files.
    """
    with zipfile.ZipFile(f'{zipname}.zip', 'w') as zipF:
        for file_name in files:
            zipF.write(
                file_name,
                compress_type=zipfile.ZIP_DEFLATED,
                )

    if remove:
        remove_files(files)

    return


def is_dimension_array(arr, ndim):
    """Assert value is array and of certain dimension."""
    valid = \
        isinstance(arr, np.ndarray) \
        and arr.ndim == ndim

    return valid


def is_higher_dimension_array(arr, ndim):
    """Assert value is array and of certain dimension."""
    valid = \
        isinstance(arr, np.ndarray) \
        and arr.ndim > ndim

    return valid


is_1d_array = partial(is_dimension_array, ndim=1)
is_2d_array = partial(is_dimension_array, ndim=2)
is_3d_array = partial(is_dimension_array, ndim=3)


def try_to_squeeze_me(arr):
    """Squeeze the arr if is array."""
    return np.squeeze(arr) if isinstance(arr, np.ndarray) else arr


def remove_files(files):
    """Remove files from disk."""
    for filename in files:
        os.remove(filename)
    return


def savetxt_w_command(fname, X, header='', fsuffix=".csv", **kwargs):
    """
    Save CSV data with info about execution command.

    Adds the command line input to the header and checks for a doubled
    defined filesuffix.
    """
    header = "{}\n{}".format(get_cli_input(), header)
    fname = "{}{}".format(fname, (not fname.endswith(fsuffix)) * fsuffix)
    np.savetxt(
        fname,
        X,
        header=header,
        fmt="%-20s",
        **kwargs)


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
