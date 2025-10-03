"""Functions for editing genomes."""

from typing import Iterable, Optional
from functools import cache
from hashlib import md5
import os

from carabiner import print_err
from joblib import Memory

from .caching import CACHE_DIR

mem = Memory(location=CACHE_DIR, verbose=0)

ATTRIBUTE_KEY_PRECEDENT = (
    "ID",
    "GeneID",
    "Name",
    "gene",
    "locus_tag",
    "old_locus_tag",
    "gene_synonym",
)


def _delete_locus(
    fasta,
    locus: Iterable[int]
):
    made_edit = False
    for seq in fasta.sequences:
        if seq.name == chr:
            seq.name = f"{seq.name}_delta-{'-'.join(interval)}"
            seq.sequence = (
                seq.sequence[:(start-1)]
                + ("N" * deletion_size)
                + seq.sequence[stop:]
            )
            made_edit = True
    if not made_edit:
        raise KeyError(f"There was no sequence called {chr} in {fasta}")
    return fasta


def _resolve_gene(
    gff,
    gene_name: str,
):
    intervals = [
        (line.columns.seqid, line.columns.start, line.columns.end, key)
        for _locus in locus
        for key in ATTRIBUTE_KEY_PRECEDENT
        for line in gff
        if key in line.attributes and line.attributes[key] == gene_name
    ]

    if len(gene_intervals) == 0:
        raise ValueError(f"Searched in {gff=}, but could not find {gene_name=}")
    print_err(f"Found {len(gene_intervals)} intervals matching {gene_name=}")
    print_err(f"Taking first match for {gene_name=}: {interval[0]}")
    return interval[0]


@cache
@mem.cache
def _resolve_gene_loci(
    gff_file: str,
    loci: Union[str, Iterable[str]]
):
    from bioino import GffFile
    if isinstance(locus, str):
        locus = [locus]
    
    gff = GffFile.from_file(gff_file)
    intervals = []
    
    for gene_name in locus:
        intervals.append(
            _resolve_gene(
                gff,
                gene_name=gene_name,
            )
        )
    intervals = sorted(intervals)
    chr, start, _ = intervals[0]
    _, _, end = intervals[-1]
    return (chr, start, end)


def delete_loci(
    fasta_file: str,
    gff_file: str,
    loci: Union[str, Union[Iterable[str], Iterable[int]]],
    cache_dir: Optional[str] = None
) -> str:

    cache_dir = cache_dir or CACHE_DIR

    if isinstance(loci, str):
        loci = []
    else:
        loci = list(loci)
    
    loci_to_delete = []
    for locus in loci:
        if (
            isinstance(locus, str) 
            or (
                isinstance(locus, (tuple, list))
                and len(locus) > 0
                and isinstance(locus[1], str)
            )
        ):  
            loci_to_delete.append(
                _resolve_gene_loci(
                    gff_file,
                    loci=locus,
                )
            )
        elif (
            isinstance(locus, (tuple, list))
            and len(locus) == 3
            and isinstance(locus[1], int)
        ):
            loci_to_delete.append(locus)
        else:
            raise ValueError(f"Invalid locus type {type(locus)} with length {len(locus)}: {locus}")

    loci_to_delete = sorted(loci_to_delete)

    _hash = md5(repr(loci_to_delete).encode()).hexdigest()
    output_file = os.path.join(cache_dir, f"{os.path.basename(fasta)}_delta-{_hash}.fna")
    
    if os.path.exists(output_file):
        return output_file
    else:
        from bioino import FastaCollection
        os.makedirs(os.path.dirname(output_file), exist_ok=True)
        fasta = FastaCollection(fasta_file)
        for locus in loci_to_delete:
            fasta = _delete_locus(fasta, locus)  
        print_err(f"Caching edited sequence at {output_file}...", end=" ")
        fasta.write(output_file)
        print_err("ok")

    return output_file
