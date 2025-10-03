"""Fetching remote data."""

from typing import Iterable, List, Optional, Union
from functools import cache
from io import BytesIO
import json
import os
from zipfile import ZipFile

from carabiner import print_err, pprint_dict

from requests import Response

from .caching import CACHE_DIR
from .http import api_get

NCBI_CACHE = os.path.join(CACHE_DIR, "ncbi")

@api_get(
    url="https://eutils.ncbi.nlm.nih.gov/entrez/eutils/espell.fcgi",
    query_key="term",
    cache_dir=NCBI_CACHE,
)
def spellcheck(query, r: Response) -> str:
    import xml.etree.ElementTree as ET
    tree = ET.parse(BytesIO(r.content))
    root = tree.getroot()
    return root.find("CorrectedQuery").text


@api_get(
    url="https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/{query}/download",
    default_params={
        "include_annotation_type": [
            "GENOME_FASTA",
            "GENOME_GFF",
        ],
        "hydrated": "FULLY_HYDRATED",
        "filename": "ncbi-dataset.zip",
    },
    cache_dir=NCBI_CACHE,
)
def download_genomic_info(
    query,
    r: Response,
    cache_dir: Optional[str] = None,
) -> List[str]:
    cache_dir = cache_dir or CACHE_DIR
    z = ZipFile(BytesIO(r.content))

    contents = z.namelist() 
    files = {
        "fasta": [
            z.extract(f, path=cache_dir) for f in contents
            if f.endswith(".fna")
        ][0],
        "gff": [
            z.extract(f, path=cache_dir) for f in contents
            if f.endswith(".gff")
        ][0],
    }

    # normalize filenames
    normalized_files = {}
    for key, f in files.items():
        _, ext = os.path.splitext(f)
        destination = os.path.join(cache_dir, f"{query}{ext}")
        print_err(f"Saving {f} at {destination}")
        os.rename(f, destination)
        normalized_files[key] = destination
    os.rmdir(os.path.dirname(f))

    return normalized_files


@api_get(
    url="https://api.ncbi.nlm.nih.gov/datasets/v2/genome/taxon/{query}/dataset_report",
    default_params={
        "filters.has_annotation": True,
        "filters.exclude_paired_reports": True,
        "filters.assembly_version": "current",
        "tax_exact_match": True,
        "table_fields": "ASSM_ACC",
    },
    cache_dir=NCBI_CACHE,
)
def taxon_to_accession(query, r: Response) -> str:
    return r.json()["reports"][0]["accession"]


@api_get(
    url="https://api.ncbi.nlm.nih.gov/datasets/v2/taxonomy/taxon_suggest/{query}",
    default_params={
        "tax_rank_filter": "species",
        "taxon_resource_filter": "TAXON_RESOURCE_FILTER_GENOME", 
    },
    cache_dir=NCBI_CACHE,
)
def name_to_taxon(query, r: Response) -> str:
    try:
        return r.json()["sci_name_and_ids"][0]["tax_id"]
    except KeyError:
        return None


def name_or_taxon_to_genome_info(
    query,
    cache_dir: Optional[str] = None,
):
    print_err(f"Fetching {query}...")
    if isinstance(query, int) or (isinstance(query, str) and query.isdigit()):
        spellchecked = str(query)
        taxon_id = spellchecked
    else:
        spellchecked = spellcheck(query)
        taxon_id = name_to_taxon(spellchecked)
    accession = taxon_to_accession(taxon_id)
    data_files = download_genomic_info(accession, cache_dir=cache_dir)
    return {
        "query": query,
        "spellchecked": spellchecked,
        "taxon_id": taxon_id,
        "accession": accession,
        "files": data_files,
    }


def fetch_landmarks(
    group: int = 0,
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
                cache_dir=cache_dir,
            ))
            pprint_dict(results[-1])
        with open(manifest_filename, "w") as f:
            json.dump(results, f, indent=4)
        
    return results


def get_landmark_ids(
    group: int = 0,
    id_keys: Optional[Iterable[Union[int, str]]] = None,
    force: bool = False,
    cache_dir: Optional[str] = None
):
    id_keys = id_keys or ("query", "taxon_id", "accession")
    landmark_info = fetch_landmarks(
        group=group,
        force=force,
        cache_dir=cache_dir,
    )
    landmarks = [
        ":".join(str(info[key]) for key in id_keys)
        for info in landmark_info
    ]
