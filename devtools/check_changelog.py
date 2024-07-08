"""
Certify the developer has input all requirements for PR.

Situations tested:

* additions are reported in CHANGELOG.rst
"""
import os
from pathlib import Path

import git


class ChangelogError(Exception):
    """Changelog error."""

    pass


folder = Path(__file__).resolve().parents[1]
changelog = Path("docs", "CHANGELOG.rst")
contributing = Path("docs", "CONTRIBUTING.rst")

repo = git.Repo(folder)

# Check if the developer is on the main branch
# We check if the head commit is the same as the main branch head commit
# as the CI will be in detached head state
if repo.active_branch.name == 'main':
    print("You are on the main branch. Nothing to check.")
    exit(0)

with open(Path(folder, changelog), "r") as fin:
    for line in fin:
        if line.startswith("v"):
            raise ChangelogError(
                "You have not updated the CHANGELOG file. "
                f"Please add a summary of your additions to {str(changelog)!r}"
                f" as described in {str(contributing)!r}."
            )
        elif line.startswith("*"):
            print("Changelog updated.")
            exit(0)
