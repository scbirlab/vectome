"""Getting and processing genome data."""

from typing import Iterable, Optional, Tuple, Union
from dataclasses import asdict, dataclass
import json
import os

from carabiner import pprint_dict, print_err

from .edits import delete_loci
from .names import _extract_species, Strain, parse_strain_label

@dataclass
class GenomeInfo:
    query: str
    spellchecked: str
    did_spellcheck: bool
    strain_info: Strain
    taxon_id: str
    accession: str
    files: Tuple[str]

    def __dict__(self):
        return asdict(self)


def name_or_taxon_to_genome_info(
    query: Union[str, int],
    check_spelling: bool = False,
    cache_dir: Optional[str] = None
):  
    from .ncbi import download_genomic_info, name_to_taxon_ncbi, spellcheck, taxon_to_accession
    print_err(f"Fetching {query}...")
    if isinstance(query, int) or (isinstance(query, str) and query.isdigit()):
        spellchecked = str(query)
        check_spelling = False
        taxon_id = spellchecked
        strain_info = parse_strain_label(taxon_id)
    else:
        species, remainder = _extract_species(query)
        spellchecked = spellcheck(species) if check_spelling else species
        strain_info = parse_strain_label(spellchecked + " " + remainder)
        search_query = strain_info.species
        for key in ("strain", "substrain"):
            if getattr(strain_info, key) is not None:
                search_query += " " + getattr(strain_info, key)
        taxon_id = name_to_taxon_ncbi(search_query, key="tax_id")
    accession = taxon_to_accession(taxon_id)
    if accession is None:
        raise KeyError(
            f"Genome lookup {taxon_id=} {search_query=} failed: {strain_info}"
        )
    data_files = download_genomic_info(
        query=accession, 
        cache_dir=cache_dir,
    )
    if (
        strain_info.deletions is not None 
        and isinstance(strain_info.deletions, list) 
        and len(strain_info.deletions) > 0
    ):
        new_fasta = delete_loci(
            fasta_file=data_files["fasta"],
            gff_file=data_files["gff"],
            loci=tuple(strain_info.deletions),
            cache_dir=cache_dir,
        )
        data_files["fasta"] = new_fasta
    return GenomeInfo(
        query=query,
        spellchecked=spellchecked,
        did_spellcheck=check_spelling,
        strain_info=strain_info,
        taxon_id=taxon_id,
        accession=accession,
        files=data_files,
    )


def fetch_landmarks(
    group: int = 0,
    check_spelling: bool = False,
    force: bool = False,
    cache_dir: Optional[str] = None
):
    from .data import load_landmarks, APPDATA_DIR

    landmarks_info = load_landmarks()

    try:
        group_queries = landmarks_info[f"group-{group}"]
    except KeyError:
        raise KeyError(
            f"Group {group} not in landmarks. Available: {', '.join(landmarks_info)}"
        )
    
    cache_dir = cache_dir or APPDATA_DIR
    cache_dir = os.path.join(cache_dir, "landmarks", f"group-{group}")
    manifest_filename = os.path.join(cache_dir, "manifest.json")

    if os.path.exists(manifest_filename) and not force:
        with open(manifest_filename, "r") as f:
            results = json.load(f)
    else:
        os.makedirs(cache_dir, exist_ok=True)

        results = []
        for q in group_queries:
            results.append(name_or_taxon_to_genome_info(
                query=q,
                check_spelling=check_spelling,
                cache_dir=cache_dir,
            ))
            pprint_dict(results[-1])
        with open(manifest_filename, "w") as f:
            json.dump(results, f, indent=4)
        
    return results


def get_landmark_ids(
    group: int = 0,
    check_spelling: bool = False,
    id_keys: Optional[Iterable[Union[int, str]]] = None,
    force: bool = False,
    cache_dir: Optional[str] = None
):
    id_keys = id_keys or ("query", "taxon_id", "accession")
    landmark_info = fetch_landmarks(
        check_spelling=check_spelling,
        group=group,
        force=force,
        cache_dir=cache_dir,
    )
    return [
        ":".join(str(info[key]) for key in id_keys)
        for info in landmark_info
    ]
