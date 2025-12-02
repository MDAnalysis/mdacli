#!/usr/bin/env python
"""Check if the CHANGELOG has been modified with respect to the main branch."""
import os

import git


class ChangelogError(Exception):
    """Changelog error."""

    pass


changelog = "docs/CHANGELOG.rst"
repo_path = os.path.realpath(os.path.join(os.path.dirname(__file__), ".."))

repo = git.Repo(repo_path)

file = repo.git.show(f"origin/main:{changelog}")

with open(os.path.join(repo_path, changelog)) as f:
    workfile = f.read()

if file.strip() == workfile.strip():
    # Check for changed files and ignore .github/ changes
    head_commit = repo.head.commit
    diff = head_commit.diff("origin/main")
    changed_files = []
    for x in diff:
        if x.a_blob.path not in changed_files:
            changed_files.append(x.a_blob.path)
        if x.b_blob is not None and x.b_blob.path not in changed_files:
            changed_files.append(x.b_blob.path)
    changed_files = [x for x in changed_files if not x.startswith(".github/")]

    if len(changed_files) > 0:
        raise ChangelogError("You have not updated the CHANGELOG file. Please "
                             f"add a summary of your additions to {changelog}.")
    else:
        print("No changes detected.")
else:
    print("CHANGELOG is up to date.")
