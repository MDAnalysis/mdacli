Contributing
============

Contributions are always welcome. Not only in the form of
code contributions via pull requests but also issues, especially
since the project is in an early development phase.

For adding new code you can follow the little introduction given below.

Fork this repository
--------------------

`Fork this repository before contributing`_.

Clone your fork
---------------

Next, clone your fork to your local machine, keep it up to date with the upstream (see later), and update the online fork with those updates.

::

    git clone https://github.com/YOUR-USERNAME/mda_cli.git
    cd mda_cli

Keep your fork up to date
-------------------------

It is important to keep your fork up-to-date with the main repository. Inside the forked project folder, every time you wish to update your fork with the main (upstream) repository do the following::

    git checkout main

    # add the upstream only the very first time
    git remote add upstream git://github.com/PicoCentauri/mda_cli.git

    git fetch upstream
    git merge upstream/main
    git pull origin main

While you are developing on a branch, you can keep it up to date with the upstream main by repeating the above commands and replacing ``main`` by the name of your development branch.

Install for developers
----------------------

If you are contributing to ``madcli``, most likely you are already a contributor of ``MDAnalysis``. Nonetheless, here are some tips.

Create a dedicated Python environment where to develop the project.

If you are using :code:`pip` follow the official instructions on `Installing packages using pip and virtual environments`_, most likely what you want is:

::

    python3 -m venv mdaclidev
    source mdaclidev/bin/activate

If you are using `Anaconda`_ go for:

::

    conda create --name mdaclidev python=3.7
    conda activate mdaclidec

Where :code:`mdaclidev` is the name you wish to give to the environment dedicated to this project.

Install ``MDAnalysis`` in your ``dev`` environment.

Either under *pip* or *conda*, install the ``mdacli`` in :code:`develop` mode, and also :ref:`tox<Uniformed Tests with tox>`.

From inside the repository main folder::

    python setup.py develop
    pip install tox

Thanks to the ``develop`` flag, any changes in the code will be automatically reflected in the installed version.

Make a new branch
-----------------

From the ``main`` branch create a new branch where to develop the new code.

::

    git checkout main
    git checkout -b new_branch


Develop the feature and keep regular pushes to your fork with comprehensible commit messages.

::

    git status
    git add (the files you want)
    git commit (add a nice commit message)
    git push origin new_branch

While you are developing, you can execute ``tox`` as needed to run unittests or inspect lint, etc. See the last section of this page.

Update CHANGELOG
----------------

Update the changelog file under :code:`docs/CHANGELOG.rst` with an explanatory bullet list of your contribution. Add that list right after the main title and before the last version subtitle::

    Changelog
    =========

    * here goes my new additions
    * explain them shortly and well

    vX.X.X (1900-01-01)
    -------------------

Also add your name to the authors list at :code:`docs/AUTHORS.rst`.

Pull Request
------------

Once you are finished, you can Pull Request you additions to the main repository, and engage with the community. Please read the ``PULLREQUEST.rst`` guidelines first, you will see them when you open a PR.

**Before submitting a Pull Request, verify your development branch passes all tests as** :ref:`described bellow<Uniformed Tests with tox>` **. If you are developing new code you should also implement new test cases.**

Also, before PR, update your development branch to the upstream main branch.

Uniformed Tests with tox
------------------------

Thanks to `Tox`_ we can have a unified testing platform where all developers are forced to follow the same rules and, above all, all tests occur in a controlled Python environment. Install ``tox`` as follows:

::

    pip install tox tox-conda
    # or
    conda install tox tox-conda -c conda-forge

You need to install ``tox-conda`` because that facilitates a lot the installation of MDAnalysis during testing.

Before creating a Pull Request from your branch, certify that all the tests pass correctly by running:

::

    tox

These are exactly the same tests that will be performed online in the Github Actions. Possibly, some tests referring to specific Python versions may fail because the interpreter is not installed. Ignored these tests.

Also, you can run individual environments if you wish to test only specific functionalities, for example:

::

    tox -e lint  # code style
    tox -e build  # packaging
    tox -e docs  # only builds the documentation
    tox -e prreq  # specific requests for PRs
    tox -e py37


.. _Tox: https://tox.readthedocs.io/en/latest/
.. _MANIFEST.in: https://github.com/MDAnalysis/mdacli/blob/main/MANIFEST.in
.. _Fork this repository before contributing: https://github.com/MDAnalysis/mdacli/network/members
.. _Pull Request: https://github.com/MDAnalysis/mdacli/pulls
.. _PULLREQUEST.rst: https://github.com/MDAnalysis/mdacli/blob/main/docs/PULLREQUEST.rst
.. _Installing packages using pip and virtual environments: https://packaging.python.org/guides/installing-using-pip-and-virtual-environments/#creating-a-virtual-environment
.. _Anaconda: https://www.anaconda.com/
