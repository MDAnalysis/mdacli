"""
Certify the developer has input all requirements for PR.

Situations tested:

* additions are reported in CHANGELOG.rst
"""
from pathlib import Path

folder = Path(__file__).resolve().parents[1]
changelog = Path('docs', 'CHANGELOG.rst')
contributing = Path('docs', 'CONTRIBUTING.rst')


class ChangelogError(Exception):
    """Changelog error."""

    pass


with open(Path(folder, changelog), 'r') as fin:
    for line in fin:
        if line.startswith('v'):
            raise ChangelogError(
                'You have not updated the CHANGELOG file. '
                f'Please add a summary of your additions to {str(changelog)!r} '
                f'as described in {str(contributing)!r}.'
                )
        elif line.startswith('*'):
            break
