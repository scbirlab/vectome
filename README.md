# 🧬 vectome

![GitHub Workflow Status (with branch)](https://img.shields.io/github/actions/workflow/status/scbirlab/vectome/python-test.yml)
![PyPI - Python Version](https://img.shields.io/pypi/pyversions/vectome)
![PyPI](https://img.shields.io/pypi/v/vectome)

**vectome** is a python package for deterministic vectorization of genomes. 

- [Installation](#installation)
- [Command-line interface](#command-line-interface)
- [Issues, problems, suggestions](#issues-problems-suggestions)
- [Documentation](#documentation)

## Installation

### The easy way

You can install the precompiled version directly using `pip`.

```bash
$ pip install vectome
```

### From source

Clone the repository, then `cd` into it. Then run:

```bash
$ pip install -e .
```

## Command-line interface

**vectome** has a command-line interface. 

```bash
$ vectome --help
```

You can generate vector embeddings by species / strain name or taxon ID.

```bash
$ vectome embed <(printf "Mycobacterium tuberculosis\n83333\nEscherichia coli CFT073")
```

The resulting vectors are based on MinHash sketches from [sourmash](), then folded into a 4096-vector using
the CountSketch method. You can make a shorter vector using e.g. `-n 1024`.

You can also deterministically project into a dense vector.

```bash
$ vectome embed <(printf "Mycobacterium tuberculosis H37Rv") --projection 16
```

Change the seed with e.g. `--seed 0`.

If you need a more interpretable vector, you can generate one based on Jaccard distances to
landmark species.

```bash
$ vectome embed <(printf "Mycobacterium tuberculosis H37Rv") --method landmark
```

Several landmark groups are available. You can set the group with `--group 0`, and 
get information about each one with `vectome info`.

```bash
$ vectome info
vectome version 0.0.1:
        group-0: {'landmarks': 113, 'manifest file': '.../vectome/vectome/data/landmarks/group-0/manifest.json', 'built': True}
        group-1: {'landmarks': 4, 'manifest file': '.../vectome/vectome/data/landmarks/group-1/manifest.json', 'built': True}
        group-2: {'landmarks': 1, 'manifest file': '.../vectome/vectome/data/landmarks/group-2/manifest.json', 'built': False}
        meta: {'cache location': '.../vectome/vectome/data/landmarks', 'cache exists': True}
```

## Issues, problems, suggestions

Add to the [issue tracker](https://www.github.com/scbirlab/vectome/issues).

## Documentation

(To come at [ReadTheDocs](https://vectome.readthedocs.org).)