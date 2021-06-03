"""Functionalities that support the creation of the command interfaces."""
import argparse
import json
from pathlib import Path


class KwargsDict(argparse.Action):
    """
    Convert input string to a dictionary.

    If string points to a ".json" file, reads the file.
    Else, attempts to convert string to dictionary using json.loads.
    """

    def __call__(self, parser, namespace, value, option_string=None):
        """Call me."""
        if value.startswith("{") and value.endswith("}"):
            try:
                jdict = json.loads(value)

            except json.decoder.JSONDecodeError as err:
                raise json.decoder.JSONDecodeError(
                    "An error ocurred when reading "
                    f"{self.dest!r} argument.",
                    err.doc,
                    err.pos,
                    ) from None
        else:
            with open(value, 'r') as fin:
                jdict = json.load(fin)

        setattr(namespace, self.dest, jdict)
