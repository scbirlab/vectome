"""Fetching remote data."""

from typing import Iterable, List, Optional, Union
from io import BytesIO
import json
import os
from zipfile import ZipFile

from carabiner import print_err

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
def name_to_taxon_ncbi(query, r: Response, key: str = "taxid", rank: Optional[str] = None) -> str:
    try:
        call_results = r.json()["sci_name_and_ids"]
    except KeyError:
        return None
    else:
        results = []
        for item in call_results:
            if rank is not None:
                try:
                    item_rank = item["rank"]
                except KeyError:
                    pass
                else:
                    if item_rank.casefold() == rank.casefold():
                        results.append(item)
            else:
                results.append(item)
        return results[0].get(key)


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
def name_to_taxon(query, r: Response, key: str = "taxonId") -> str:
    try:
        return r.json()["results"][0]["taxonomy"][key]
    except KeyError:
        return None
