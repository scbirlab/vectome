"""Getting and processing genome data."""

from typing import Iterable, Optional, Tuple, Union
from dataclasses import asdict, dataclass
from functools import cache, partial
import json
import os

from carabiner import pprint_dict, print_err

from .data.download import download_landmark_cache
from .edits import delete_loci
from .names import _extract_species, Strain, parse_strain_label
from . import __version__

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

    @classmethod
    def from_name_or_taxon(
        cls,
        query: Union[str, int],
        check_spelling: bool = False,
        quiet: bool = False,
        cache_dir: Optional[str] = None,
        strict: bool = False,    # `True` allows fallback to species name
        _landmark: bool = False  # prevents cache hits on landmark downloads
    ) -> 'GenomeInfo':
        return name_or_taxon_to_genome_info(
            query=query,
            check_spelling=check_spelling,
            cache_dir=cache_dir,
            strict=strict,
            quiet=quiet,
            _landmark=_landmark,
        )

    @classmethod
    def from_file(
        cls,
        file: str,
        gff: Optional[str] = None,
        quiet: bool = False,
        cache_dir: Optional[str] = None,
        _landmark: bool = False  # prevents cache hits on landmark downloads
    ) -> 'GenomeInfo':
        return file_to_genome_info(
            file=file,
            gff=gff,
            quiet=quiet,
            cache_dir=cache_dir,
            _landmark=_landmark,
        )

    @classmethod
    def from_any(
        cls,
        query: Union[str, int],
        check_spelling: bool = False,
        cache_dir: Optional[str] = None,
        strict: bool = False,    # `True` allows fallback to species name
        quiet: bool = False,
        _landmark: bool = False  # prevents cache hits on landmark downloads
    ) -> 'GenomeInfo':
        if str(query).startswith(("file://", "url://", "http://", "https://")):
            return cls.from_file(
                file=query,
                gff=None,
                quiet=quiet,
                cache_dir=cache_dir,
                _landmark=_landmark,
            )
        else:
            return cls.from_name_or_taxon(
                query=query,
                check_spelling=check_spelling,
                quiet=quiet,
                cache_dir=cache_dir,
                strict=strict,
                _landmark=_landmark,
            )


def resolve_file_url(
    file: str, 
    cache_dir: str, 
    quiet: bool = False
) -> Tuple[bool, str]:
    if file.startswith("file://"):
        file = file.split("file://")[-1]
        is_url = False
    elif file.startswith("url://"):
        file = "https://" + file.split("url://")[-1]
        is_url = True
    elif file.startswith(("http://", "https://")):
        is_url = True
    else:
        is_url = False
    if is_url:
        from .http import download_url
        dl_dir = os.path.join(cache_dir, "downloads")
        os.makedirs(dl_dir, exist_ok=True)
        destination = os.path.join(dl_dir, os.path.basename(file))
        
        return is_url, download_url(
            file, 
            destination=destination,
            cache_dir=cache_dir,
            quiet=quiet,
        )
    elif os.path.exists(file):
        return is_url, file
    else:
        raise FileNotFoundError(f"File does not exist or is not URL: {file}")


@cache
def file_to_genome_info(
    file: str,
    gff: Optional[str] = None,
    quiet: bool = False,
    cache_dir: Optional[str] = None,
    _landmark: bool = False  # prevents cache hits on landmark downloads
):  
    from .ncbi import _normalize_and_compress
    
    is_url, fasta_file = resolve_file_url(file, cache_dir=cache_dir)
    data_files = {
        "fasta": fasta_file,
        "gff": resolve_file_url(gff, cache_dir=cache_dir) if gff is not None else None,
    }
    if is_url:
        data_files = _normalize_and_compress(
            data_files,
            query=os.path.splitext(os.path.basename(file))[0],
            cache_dir=cache_dir,
            move=is_url,
            quiet=quiet,
        )
    return GenomeInfo(
        query=file,
        spellchecked=file,
        did_spellcheck=False,
        strain_info=None,
        taxon_id=None,
        accession=None,
        files=data_files,
    )


@cache
def name_or_taxon_to_genome_info(
    query: Union[str, int],
    check_spelling: bool = False,
    quiet: bool = False,
    hide_progress: bool = False,
    cache_dir: Optional[str] = None,
    strict: bool = False,    # `True` allows fallback to species name
    _landmark: bool = False  # prevents cache hits on landmark downloads
):  
    from .ncbi import download_genomic_info, name_to_taxon_ncbi, spellcheck, taxon_to_accession
    if query is None:
        return None
    if not quiet:
        print_err(f"Fetching {query}...")
    if isinstance(query, int) or (isinstance(query, str) and query.isdigit()):
        spellchecked = str(query)
        check_spelling = False
        taxon_id = spellchecked
        strain_info = parse_strain_label(taxon_id)
        search_query = strain_info.species
    elif isinstance(query, str) and query.startswith(("GCF_", "GCA_")) and query[-1].isdigit():
        spellchecked = query
        check_spelling = False
        taxon_id = None
        strain_info = parse_strain_label(query)
        accession = query
    else:
        species, remainder = _extract_species(query)
        spellchecked = spellcheck(species) if check_spelling else str(species)
        strain_info = parse_strain_label(spellchecked + " " + remainder)
        search_query = strain_info.species
        for key in ("strain", "substrain"):
            if getattr(strain_info, key) is not None:
                search_query += " " + getattr(strain_info, key)
        taxon_id = name_to_taxon_ncbi(search_query, key="tax_id")
    accession = taxon_to_accession(taxon_id)
    if accession is None:
        if strict or _landmark:
            raise KeyError(
                f"Genome lookup {taxon_id=} {search_query=} failed: {strain_info}"
            )
        else:
            return name_or_taxon_to_genome_info(
                query=strain_info.species,
                check_spelling=check_spelling,
                cache_dir=cache_dir,
                strict=True,
                _landmark=_landmark,
            )
    if not quiet:
        print_err(f"[INFO] Parsed {search_query=} -> {taxon_id=}")
    data_files = download_genomic_info(
        query=accession,
        quiet=quiet,
        cache_dir=cache_dir,
        _landmark=_landmark,
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
    quiet: bool = False,
    hide_progress: bool = False,
    allow_missing_files: bool = False,
    redownload: bool = False,
    cache_dir: Optional[str] = None
):
    from tqdm.auto import tqdm
    from requests import HTTPError

    from .data import load_landmarks, APPDATA_DIR

    landmarks_info = load_landmarks()

    try:
        group_queries = landmarks_info[f"group-{group}"]
    except KeyError:
        raise KeyError(
            f"Group {group} not in landmarks. Available: {', '.join(landmarks_info)}"
        )
    
    cache_dir = cache_dir or APPDATA_DIR
    landmark_version = os.environ.get('VECTOME_LANDMARKS_VERSION', f"v{__version__}")
    landmarks_dir = os.path.join(cache_dir, "landmarks", landmark_version, f"group-{group}")
    manifest_filename = os.path.join(landmarks_dir, "manifest.json")

    if os.path.exists(manifest_filename) and not force:
        with open(manifest_filename, "r") as f:
            try:
                results = json.load(f)
            except json.JSONDecodeError as e:
                print_err(f.read())
                return fetch_landmarks(
                    group=group,
                    check_spelling=check_spelling,
                    force=True,
                    quiet=quiet,
                    hide_progress=hide_progress,
                    cache_dir=cache_dir,
                )
    elif not redownload:
        try:
            _ = download_landmark_cache(
                suffix=group,
                quiet=quiet,
                cache_dir=cache_dir,
            )
        except HTTPError:
            print_err("Could not find downloadable cache. Downloading from original source")
            return fetch_landmarks(
                group=group,
                check_spelling=check_spelling,
                force=force,
                redownload=True,
                allow_missing_files=allow_missing_files,
                quiet=quiet,
                hide_progress=hide_progress,
                cache_dir=cache_dir,
            )
        else:
            with open(manifest_filename, "r") as f:
                results = json.load(f)
    else:
        os.makedirs(landmarks_dir, exist_ok=True)

        results = []
        errors = {}
        _iter = iter if hide_progress else partial(tqdm, desc="Fetching landmarks") 
        for q in _iter(group_queries):
            try:
                genome_info = GenomeInfo.from_any(
                    query=q,
                    check_spelling=check_spelling,
                    cache_dir=landmarks_dir,
                    _landmark=True,
                )
            except Exception as e:
                genome_info = None
                errors[q] = e
                if not quiet:
                    print_err(e)
                    print_err(f"[WARN] Failed to get genome info for query {q}!")
                # raise e
            else:
                if not quiet:
                    pprint_dict(asdict(genome_info), message="Parsed strain name:")
            results.append(asdict(genome_info))
        if len(errors) > 0:
            message = f"[ERROR] Failed to fetch {len(errors)} queries!"
            print_err(message)
            print_err("\n".join(errors))
            raise ValueError(errors[list(errors)[0]])
        with open(manifest_filename, "w") as f:
            json.dump(results, f, indent=4)

    # check all files exist, otherwise delete manifest and regenerate
    rebuild = False
    for item in results:
        for key, filename in item["files"].items():
            if not allow_missing_files and not os.path.exists(filename):
                if not quiet:
                    print_err(
                        f"[WARN] The '{key}' file ({filename}) for {item['query']} is missing!",
                        f"Deleting manifest and rebuilding group {group} landmarks.",
                    )
                os.remove(manifest_filename)
                rebuild = True
                break
            else:
                if not quiet:
                    print_err(
                        f"[INFO] Found '{key}' file ({filename}) for {item['query']}",
                    )

    if rebuild:
        return fetch_landmarks(
            group=group,
            check_spelling=check_spelling,
            force=True,
            redownload=redownload,
            allow_missing_files=allow_missing_files,
            quiet=quiet,
            hide_progress=hide_progress,
            cache_dir=cache_dir,
        )
    else:
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
        ":".join(
            os.path.basename(str(info[key])) 
            for key in id_keys if key is not None
        )
        for info in landmark_info
    ]
