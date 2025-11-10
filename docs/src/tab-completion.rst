===================
Tab-Completion
===================

``mdacli`` includes built-in support for command-line tab-completion using

Activation
==========

The activation method depends on your shell:

Bash
----

**Temporary (current session only)**::

    eval "$(register-python-argcomplete mda)"

**Permanent (recommended)**

Add to your ``~/.bashrc``::

    echo 'eval "$(register-python-argcomplete mda)"' >> ~/.bashrc
    source ~/.bashrc

Zsh
---

Add to your ``~/.zshrc``::

    autoload -U bashcompinit
    bashcompinit
    eval "$(register-python-argcomplete mda)"

Then reload::

    source ~/.zshrc

Fish
----

Generate the completion file::

    register-python-argcomplete --shell fish mda > ~/.config/fish/completions/mda.fish

Restart your Fish shell or run::

    source ~/.config/fish/config.fish

Tcsh
----

Add to your shell startup file::

    eval `register-python-argcomplete --shell tcsh mda`

Usage Examples
==============

Once enabled, tab-completion works for:

**Module names**::

    mda <TAB>
    # Shows: AlignTraj, AverageStructure, Contacts, DensityAnalysis, ...

**Partial module names**::

    mda RM<TAB>
    # Shows: RMSD, RMSF

**Options and flags**::

    mda RMSD -<TAB>
    # Shows: -s, -f, -atomgroup, -b, -e, -dt, -v, --debug, --version, ...

**Case insensitive**::

    mda rmsd<TAB>    # Also works
    mda RmSd<TAB>    # Also works

Troubleshooting
===============

Tab-completion not working
--------------------------

1. **Verify argcomplete is installed**::

       python -c "import argcomplete; print(argcomplete.__version__)"

2. **Check if activation command was added**::

       grep "register-python-argcomplete mda" ~/.bashrc

3. **Reload your shell**::

       source ~/.bashrc  # or restart terminal

4. **Test basic completion**::

       mda <TAB>

Still not working
-----------------

- Make sure you've restarted your terminal or sourced the configuration file
- For Zsh, ensure ``bashcompinit`` is loaded before argcomplete
- Check that ``mda`` is in your PATH: ``which mda``
- Try running the registration command manually in your current shell

Global activation (for all Python scripts)
-------------------------------------------

To enable argcomplete for all Python scripts at once::

    activate-global-python-argcomplete

This requires root/admin privileges and will enable completion system-wide.
