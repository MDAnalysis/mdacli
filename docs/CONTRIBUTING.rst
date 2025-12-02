Contributing
============

Contributing to this project is easy because we have set up a powerful
CI environment. :relaxed: This document will guide you step-by-step.

Contributions are always welcome. Not only in the form of
code contributions via pull requests but also issues, especially
since the project is in an early development phase.

The way to interact with us is to `fork` the `mdacli` repository to your
account, create a development branch in your fork, and finally pull
request your changes to the main repository for review. Next, we present
you guidelines for this process. If you are a `git` pro you may wish to
apply your own `git` method; if not, you are safe following ours.

Fork this repository
-------------------------------------------------------------------------

`Fork this repository before contributing`_.

Clone the main project to your computer
-------------------------------------------------------------------------

Next, clone the main repository to your local machine.::

    git clone https://github.com/MDAnalysis/mdacli.git
    cd mdacli

Add your fork as remote
-------------------------------------------------------------------------

You can't create a branch in the main repository :no_entry_sign:, you
need to make it in your fork :smile: , first you need to add your fork
as a `remote` in the cloned repository (following the previous commands)
::

    git checkout main
    git remote add myfork git://github.com/<your-username>/mdacli.git
    git fetch myfork


Create a new branch and start developing
-------------------------------------------------------------------------

Now create a new branch and start developing::

    git checkout -b <my-new-branch-with-a-nice-name>

Develop your code. You should commit your changes to your fork in
encapsulated steps. That is, when you finish doing some task, you should
commit it.

    git add <the new files>
    git commit -m "<your good commit message>"
    git push myfork <my-new-branch-with-a-nice-name>

You will see that your changes are now in your fork and branch.

Pull Request your changes
----------------------------------

Once you finish your changes create a Pull Request to the main `mdacli`
repository so we can review your contribution, give you feedback, and
hopefully accept it :relaxed:. However, before PR, continue reading this
guideline!

Install for developers
----------------------

If you are contributing to ``madcli``, most likely you are already a
contributor of ``MDAnalysis v2``. If you already have `MDAnalysis v2`
installed, install `mdacli` in the same Python environment. For that,
inside `mdacli`'s repository::

    python setup.py develop --no-deps

It is very important to use the `develop` flag, so that the changes you
make in the code are always reflected (real time) in your installation.
The `--no-deps` avoids installing `MDAnalysis` twice.

If you don't have `MDAnalysis v2` installed how did you run `mdacli` in
the first place? Please refer to :ref:`Installation` before continuing.

The whole `mdacli` Continuous Integration pipeline is based on `Tox_`.
So you need to install `tox` to be our friend :smile_cat: ::

    pip install tox tox-conda
    # or
    conda install tox tox-conda -c conda-forge


Running tests with tox
---------------------------

Before creating a Pull Request from your branch, certify that all the
tests pass correctly by running:

::

    tox

These are exactly the same tests that will be performed online in our
Github Actions workflows. Possibly, some tests referring to specific
Python versions may fail because the interpreter is not installed;
ignored these tests.

Also, you can run individual environments if you wish to test only
specific functionalities, for example:

::

    tox -e lint  # code style
    tox -e build  # packaging
    tox -e docs  # only builds the documentation
    tox -e test  # testing


Update CHANGELOG
----------------

Update the changelog file under :code:`docs/CHANGELOG.rst` with an
explanatory bullet list of your contribution bellow the `CHANGELOG`
title. Add that list right after the main title and before the last
version subtitle::

    Changelog
    =========

    * here goes my new additions
    * explain them shortly and well

    vX.X.X (1900-01-01)
    -------------------

Also add your name to the authors list at :code:`docs/AUTHORS.rst`.

Pull Request
------------

Once you are finished, you can Pull Request you additions to the main
repository and engage with the community. Please read the `docs/PULLREQUEST.rst`
guidelines first, you will see them when you open a PR.

**Before submitting a Pull Request, verify your development branch
passes all tests as** :ref:`described<Running tests with tox>` **. If
you are developing new code you should also implement new test cases.**

Also, before PR, update your development branch to the upstream main
branch to certify there are no incompatibilities::

    git checkout main
    git pull
    git checkout <my-new-branch-with-a-nice-name>
    git merge --no-ff main


Correct any conflicts that may appear. It there are no conflicts, you
are good to go (Pull Request).

.. _Tox: https://tox.readthedocs.io/en/latest/
.. _Fork this repository before contributing: https://github.com/MDAnalysis/mdacli/network/members

How to Perform a Release
-------------------------
1. **Prepare a Release Pull Request**

   - Based on the main branch create branch ``release-0.1.35`` and a PR.
   - Ensure that all `CI tests <https://github.com/MDAnalysis/mdacli/actions>`_ pass.
   - Optionally, run the tests locally to double-check.

2. **Update the Changelog**

   - Edit the changelog located in ``docs/CHANGELOG.rst``:
      - Add a new section for the new version, summarizing the changes based on the
        PRs merged since the last release.
      - Leave a placeholder section titled *Unreleased* for future updates.

3. **Merge the PR and Create a Tag**

   - Merge the release PR.
   - Update the ``main`` branch and check that the latest commit is the release PR with
     ``git log``
   - Create a tag on directly the ``main`` branch.
   - Push the tag to GitHub. For example for a release of version ``2025.3``:

     .. code-block:: bash

        git checkout main
        git pull
        git tag -a v0.1.35 -m "Release v0.1.35"
        git push --tags

4. **Finalize the GitHub Release**

   - Once the PR is merged, the CI will automatically:
      - Publish the package to PyPI.
      - Create a draft release on GitHub.
   - Update the GitHub release notes by pasting the changelog for the version.

5. **Merge Conda Recipe Changes**

   - May resolve and then merge an automatically created PR on the `conda recipe
     <https://github.com/conda-forge/mdacli-feedstock>`_.
   - Once thus PR is merged and the new version will be published automatically on the
     `conda-forge <https://anaconda.org/conda-forge/mdacli>`_ channel.
