"""Functions for editing genomes."""

from typing import Iterable, Optional, Union
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
    chr, start, stop = locus
    if start == stop and start < 0:
        return fasta  # interval wasn't found
    made_edit = False
    deletion_size = stop - start + 1
    for seq in fasta.sequences:
        try:
            this_chr = seq.name
        except AttributeError as e:
            raise AttributeError(f"Sequence is type {type(seq)}: {seq}")
        if this_chr == chr:
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
    strict: bool = False
):
    intervals = [
        (line.columns.seqid, line.columns.start, line.columns.end, key)
        for key in ATTRIBUTE_KEY_PRECEDENT
        for line in gff.lines
        if key in line.attributes and line.attributes[key].casefold() == gene_name.casefold()
    ]

    if len(intervals) == 0:
        if strict:
            raise ValueError(f"Searched in GFF, but could not find {gene_name=}")
        else:
            print_err(f"[WARN] Could not find feature {gene_name=} in GFF")
            line = list(gff.lines)[0]
            return (line.columns.seqid, -1, -1, ATTRIBUTE_KEY_PRECEDENT[0])
    print_err(f"Found {len(intervals)} intervals matching {gene_name=}")
    interval = intervals[0]
    print_err(f"Taking first match for {gene_name=}: {interval}")
    return interval


@cache
@mem.cache
def _resolve_gene_loci(
    gff_file: str,
    loci: Union[str, Iterable[str]]
):
    from bioino import GffFile
    if isinstance(loci, str):
        loci = [loci]
    
    gff = GffFile.from_file(gff_file)
    gff.lines = tuple(gff.lines)
    intervals = []
    
    for gene_name in loci:
        intervals.append(
            _resolve_gene(
                gff,
                gene_name=gene_name,
            )
        )
    intervals = sorted(intervals)
    chr, start, _ = intervals[0][:3]
    _, _, end = intervals[-1][:3]
    return (chr, start, end)


def delete_loci(
    fasta_file: str,
    gff_file: str,
    loci: Union[str, Union[Iterable[str], Iterable[int]]],
    cache_dir: Optional[str] = None
) -> str:

    cache_dir = cache_dir or CACHE_DIR

    if isinstance(loci, str):
        loci = [loci]
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
            and len(locus) >= 3
            and isinstance(locus[1], int)
        ):
            loci_to_delete.append(locus)
        else:
            raise ValueError(f"Invalid locus type {type(locus)} with length {len(locus)}: {locus}")

    loci_to_delete = sorted(loci_to_delete)

    fasta_basename = os.path.basename(fasta_file)
    _hash = md5((repr(loci_to_delete) + fasta_basename).encode()).hexdigest()
    output_file = os.path.join(cache_dir, f"{fasta_basename.split('_delta-')[0]}_delta-{_hash}.fna")
    
    if os.path.exists(output_file):
        return output_file
    else:
        from bioino import FastaCollection
        os.makedirs(os.path.dirname(output_file), exist_ok=True)
        fasta = FastaCollection.from_file(fasta_file)
        fasta.sequences = tuple(fasta.sequences)
        for locus in loci_to_delete:
            print_err(f"Deleting {locus}...")
            fasta = _delete_locus(
                fasta=fasta, 
                locus=locus,
            )  
        print_err(f"Caching edited sequence at {output_file}...", end=" ")
        with open(output_file, "w") as fh:
            fasta.write(fh)
        print_err("ok")

    return output_file
