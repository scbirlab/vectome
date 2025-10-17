"""Fetching remote data."""

from typing import Iterable, List, Optional, Union
from io import BytesIO
import json
import os

from carabiner import print_err

from .caching import CACHE_DIR
from .http import api_get

NCBI_CACHE = os.path.join(CACHE_DIR, "ncbi")


@api_get(
    url="https://eutils.ncbi.nlm.nih.gov/entrez/eutils/espell.fcgi",
    query_key="term",
    cache_dir=NCBI_CACHE,
)
def spellcheck(query, r) -> str:
    import xml.etree.ElementTree as ET
    tree = ET.parse(BytesIO(r.content))
    root = tree.getroot()
    try:
        return root.find("CorrectedQuery").text
    except Exception:
        return None


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
    r,
    cache_dir: Optional[str] = None,
    _landmark: bool = False  # prevents cache hits on landmark downloads
) -> List[str]:

    from zipfile import ZipFile
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
    if all(os.path.exists(f) for key, f in normalized_files.items()):
        return normalized_files
    else:
        raise IOError(f"Some files are missing! {({key: f for key, f in normalized_files.items() if not os.path.exists(f)})}")


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
def taxon_to_accession(query, r) -> str:
    call_results = r.json().get("reports")
    if call_results is not None and isinstance(call_results, list) and len(call_results) > 0:
        return call_results[0].get("accession")
    else:
        return None


@api_get(
    url="https://api.ncbi.nlm.nih.gov/datasets/v2/taxonomy/taxon_suggest/{query}",
    default_params={
        "tax_rank_filter": "species",
        "taxon_resource_filter": "TAXON_RESOURCE_FILTER_GENOME", 
    },
    cache_dir=NCBI_CACHE,
)
def name_to_taxon_ncbi(query, r, key: str = "tax_id", rank: Optional[str] = None) -> str:
    call_results = r.json().get("sci_name_and_ids")
    if call_results is not None and isinstance(call_results, list):
        if rank is None:
            rank_results = call_results
        else:
            rank_results = []
            for item in call_results:
                try:
                    item_rank = item["rank"]
                except KeyError:
                    pass
                else:
                    if isinstance(item_rank, str) and item_rank.casefold() == rank.casefold():
                        rank_results.append(item)
        if len(rank_results) > 0:
            return rank_results[0].get(key)
    
    return None


@api_get(
    url="https://rest.uniprot.org/proteomes/search",
    default_params={
        "size": 1,
        "fields": ["organism", "organism_id"],
        "sort": "organism_name asc",
    },
    query_key="query",
    cache_dir=NCBI_CACHE,
)
def name_to_taxon(query, r, key: str = "taxonId") -> str:
    call_results = r.json().get("results")
    if call_results is not None and isinstance(call_results, list) and len(call_results) > 0:
        return call_results[0].get("taxonomy")
    else:
        return None
