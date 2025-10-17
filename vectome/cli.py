"""Command-line interface for vectome."""

from argparse import FileType, Namespace
from functools import cache
import sys

from carabiner import (
    pprint_dict,
    print_err,
    clicommand,
    CLIApp,
    CLICommand,
    CLIOption,
)

from . import app_name, __author__, __version__
from .caching import CACHE_DIR


@clicommand(message="Running embed")
def embed(args: Namespace) -> None:
    from .vectorize import vectorize
    strains = tuple(line.rstrip() for line in args.strain)
    vectorize = cache(vectorize)
    vectors = vectorize(
        strains,
        check_spelling=args.spellcheck,
        method=args.method,
        group=args.group,
        dim=args.dimensionality,
        projection=args.projection,
        seed=args.seed,
        cache_dir=args.cache,
    )
    
    if args.method == "landmark" and args.projection is None:
        from .genomes import get_landmark_ids
        header = get_landmark_ids(
            group=args.group, 
            cache_dir=args.cache,
        )
    else:
        header = list(map(str, range(vectors.shape[1])))
    
    if header is not None:
        print(
            "\t".join(["query"] + header),
            file=args.output,
        )
    for row, query in zip(vectors, strains):
        print(
            "\t".join([query] + list(map(str, row))), 
            file=args.output,
        )
    return None


@clicommand(message="Running build")
def build(args: Namespace) -> None:
    from .sketching import sketch_landmarks
    results = sketch_landmarks(
        group=args.group,
        check_spelling=args.spellcheck,
        force=args.force,
        cache_dir=args.cache,
    )
    return None


@clicommand(message="Running info")
def info(args: Namespace) -> None:
    from .data import landmark_info
    pprint_dict(
        landmark_info(cache_dir=args.cache),
        message=f"vectome version {__version__}"
    )
    return None


def main() -> None:

    """Main function for CLI app.

    """

    cache_dir = CACHE_DIR

    options = {
        "group": CLIOption(
            "group",
            type=int,
            help="Landmark group number.",
        ),
        "group flag": CLIOption(
            "--group", "-g",
            type=int,
            default=0,
            help="Landmark group number.",
        ),
        "strain": CLIOption(
            "strain",
            type=FileType("r"),
            nargs='?',
            default=sys.stdin,
            help="File listing name or accession ID to vectorize, one per line. Default: STDIN",
        ),
        "dimensionality": CLIOption(
            "--dimensionality", "-n",
            type=int,
            default=None,
            help="Vector size for hash folding. Default: 4096",
        ),
        "method": CLIOption(
            "--method", "-m",
            type=str,
            default="countsketch",
            choices=["countsketch", "landmark"],
            help="Vectorizaton method.",
        ),
        "force": CLIOption(
            "--force", "-f",
            action="store_true",
            help="Ignore cache and rebuild.",
        ),
        "spellcheck": CLIOption(
            "--spellcheck", "-s",
            action="store_true",
            help="Ignore cache and rebuild.",
        ),
        "projection": CLIOption(
            "--projection", "-z",
            type=int,
            default=None,
            help="If provided, the vector will be projected down to this dimensionality. Default: no projection.",
        ),
        "seed": CLIOption(
            "--seed", "-i",
            type=int,
            default=42,
            help="Seed for reproducible randomness.",
        ),
        "cache": CLIOption(
            "--cache", "-c",
            type=str,
            default=cache_dir,
            help="Where to cache data.",
        ),
        "output": CLIOption(
            "--output", "-o",
            default=sys.stdout,
            type=FileType("w"),
            help="Output file. Default: STDOUT",
        ),
    }

    app = CLIApp(
        app_name,
        version=__version__,
        description="Deterministic vectorization of genomes.",
        commands=[
            CLICommand(
                "embed",
                main=embed,
                description="Turn one or more species/strain names or accession IDs into a vector.",
                options=[
                    options["strain"],
                    options["spellcheck"],
                    options["dimensionality"],
                    options["method"],
                    options["projection"],
                    options["seed"],
                    options["group flag"],
                    options["cache"],
                    options["output"],
                ],
            ),
            CLICommand(
                "build",
                main=build,
                description="Build landmark database.",
                options=[
                    options["group"],
                    options["spellcheck"],
                    options["force"],
                    options["cache"].replace(default=None),
                ],
            ),
            CLICommand(
                "info",
                main=info,
                description="Provide information on landmarks.",
                options=[
                    options["cache"].replace(default=None),
                ],
            ),
        ],
    )
    app.run()
    return None


if __name__ == '__main__':
    main()
