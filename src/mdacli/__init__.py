"""Command line interfaces for MDAnalysis analysis classes."""

from .cli import maincli, setup_clients

__version__ = '0.0.0'


def main():
    """Execute main CLI entry point."""
    maincli(setup_clients())
