#!/usr/bin/env python3
# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
#
# Copyright (c) 2021 Authors and contributors
#
# Released under the GNU Public Licence, v2 or any higher version
# SPDX-License-Identifier: GPL-2.0-or-later

"""Useful helper functions for running the cli."""

import inspect
import re
import sys
from collections import defaultdict

import MDAnalysis as mda


def _exit_if_a_is_b(obj1, obj2, msg):
    """Exit if `obj1` and `obj2` are the same."""
    if obj1 is obj2:
        sys.exit(msg)


def split_time_unit(s):
    """
    Split time and units.

    Follows the regex:: https://regex101.com/r/LZAbil/2

    Returns
    -------
    tuple (float, str)
        Value as float, units as str.

    Raises
    ------
    IndexError
        Tuple could not be found. This happens when a number is not
        present in the start of the string.
    """
    type_regex = re.compile(r'^(\-?\d+\.?\d*|\-?\.\d+|\-?\.?\d+[eE]\-?\d+|-?\d+\.?\d*[eE]\d+)($|[a-z]*$)')  # noqa: E501
    value, unit = type_regex.findall(s)[0]
    return float(value), unit


def convert_str_time(x, dt, base_unit='ps'):
    """
    Convert a string `x` into a frame number based on given `dt`.

    If `x` does not contain any units its assumed to be a frame number
    already.

    See :func:`split_time_unit`.

    Parameters
    ----------
    x : str
        the input string
    dt : float
        the time step in ps

    Returns
    -------
    int
        frame number

    Raises
    ------
    ValueError
        The input does not contain any units but is not an integer.
    """
    val, unit = split_time_unit(x)
    if unit != "":
        val = mda.units.convert(val, unit, "ps")
        return int(val // dt)
    elif val % 1 != 0:  # the number is not int'able
        raise ValueError(
            "Only integers or time step combinations (´12ps´) "
            "are valid for frame selection"
            )
    else:
        return int(val)


def parse_callable_signature(callable_obj):
    """
    Parse a callable signature to a convenient dictionary for CLI creation.

    The parameters used in the CLI are a combination of the callable
    signature and the information in the callable docstring.

    Parameters
    ----------
    callable_obj : callable
        The callable object to inspect. Details of this object required
        for the creation of a CLI are added to the `storage_dict`.

    Returns
    -------
    dict
    """
    storage_dict = {}
    storage_dict["callable"] = callable_obj

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
                # else reaches if the parameter in the signature is not present
                # in the docstring. It shouldn't, but just in case :-)
                # unless we explicitly decide not to consider any parameters
                # not referenced in the documentation, this should be kept
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
                        # type taken form docstring
                        "type": doc_param['type'],
                        # but default taken from signature
                        "default": sig_param.default,
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
    storage_dict["positional"] = positional_args
    storage_dict["optional"] = optional_args
    storage_dict["desc"] = summary
    storage_dict["desc_long"] = summary_extended

    return storage_dict


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
    type_regex = re.compile(r'^(\w+|\{.*\})|(?<=\`\~)(.*?)(?=\`)')

    # goes back to front to register descriptions first ;-)
    # considers only the Parameters section
    for line in doc_lines[par_i: end_param_line][::-1]:
        if ' : ' in line:
            par_name, others_ = line.split(' : ')
            par_type = [_ for _ in type_regex.findall(others_)[0] if _][0]
            params[par_name]['type'] = par_type
            params[par_name]['desc'] = ' '.join(desc_tmp[::-1])
            desc_tmp.clear()
        else:
            desc_tmp.append(line)

    return summary, summary_extended, params
