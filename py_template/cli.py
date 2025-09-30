"""Command-line interface for py_template."""

from . import __version__


def main() -> None:

    """Main function for CLI app.

    Put argument parsing, etc here.

    Examples
    ========
    >>> 1 + 1
    2

    """

    import sys
    print(f"Running app version {__version__}", file=sys.stderr)
    return None


if __name__ == '__main__':
    main()
