from typing import List, Optional, Tuple, Union
from dataclasses import dataclass, field
import re

from .ncbi import name_to_taxon_ncbi

# --- Helpers ---
@dataclass
class Strain:
    query: str
    species: str
    remainder: str
    strain: Optional[str] = field(default=None)
    substrain: Optional[str] = field(default=None)
    deletions: Optional[List[Union[str, List[str]]]] = field(default_factory=list)
    mutations: Optional[List[str]] = field(default_factory=list)


def _strip_punctuation(s: str) -> str:
    return s.strip(" ,;:/()[]{}")


def _split_operon(query: str) -> List[str]:
    """Split operon shorthand like 'acrAB' into ['acrA','acrB'].

    Returns
    =======
    list[str]
        Gene names.

    """
    m = re.fullmatch(r"([a-z]{3})([A-Z]{2,})", query)
    if not m:
        return [query]
    base, caps = m.groups()
    return [base + c for c in caps]


def _extract_species(
    query: str,
    normalize: bool = False
) -> Tuple[str, str]:
    """Extract 'Genus species' (with optional abbreviated genus) from the start.

    Returns
    =======
    tuple[str]
        (species_full, remainder_text)

    """
    query = query.strip()
    # e.g., "E. coli", "Escherichia coli", allow extra words after species
    m = re.match(
        r"^\s*(?P<genus>[A-Z][a-z]*|[A-Z]\.)\s+(?P<species>[a-z][a-z_-]+)\b", 
        query,
    )
    if not m:
        return query, ""
    genus = m.group("genus")
    species = m.group("species").replace("_", "-")
    remainder = query[m.end():].strip()
    species_full = f"{genus} {species}"
    if normalize:
        species_full = name_to_taxon_ncbi(species_full, key="sci_name", rank="species")
        species_full = _extract_species(species_full, normalize=False)[0]
    return species_full, remainder


def _extract_strain_and_substrain(query: str) -> Tuple[Optional[str], Optional[str], str]:
    """
    Try to pull 'strain' and 'substrain' from phrases like:
    'K-12 substr. MG1655', 'subsp. H37Rv', 'strain ATCC 12345', 'K-12 MG1655'

    Returns 
    =======
    tuple
        (strain, substrain, remainder_text).
    """
    strain, substrain = None, None

    # substrain indicators
    m = re.search(
        r"\b(?:substr(?:ain)?\.?|subsp\.?|variant)\s+([A-Za-z0-9._-]+)", 
        query, 
        flags=re.IGNORECASE,
    )
    if m:
        substrain = _strip_punctuation(m.group(1))
        query = (query[:m.start()] + query[m.end():]).strip()

    # strain indicators
    m = re.search(
        r"\b(?:strain|str\.|serovar\.?)\s+(ATCC\s[0-9]+|[A-Za-z0-9._-]+)\b", 
        query, 
        flags=re.IGNORECASE,
    )
    if m:
        strain = _strip_punctuation(m.group(1))
        query = (query[:m.start()] + query[m.end():]).strip()

    # fallback: common alphanumeric early in the string
    if strain is None:
        m = re.search(
            r"\b([A-Za-z]\-[0-9]+|[A-Za-z0-9]{2,})\b", 
            query,
        )
        if m:
            candidate = _strip_punctuation(m.group(1))
            # Heuristic: treat as strain if it looks like K-12 or a lab nickname, and is not a gene token
            if re.match(r"^[A-Za-z]\-\d+$", candidate) or re.match(r"^[A-Z][A-Za-z0-9._-]{2,}$", candidate):
                strain = candidate
                query = (query[:m.start()] + query[m.end():]).strip()

    # second token as substrain if we saw two alphanum tokens in a row (e.g., "K-12 MG1655")
    if substrain is None:
        m = re.search(r"\b([A-Za-z0-9._-]{3,})\b", query)
        if m:
            candidate = _strip_punctuation(m.group(1))
            if re.match(r"^[A-Z]{0,2}\d*[A-Za-z0-9._-]+$", candidate) and not re.match(r"^[a-z]+[A-Z]{2,}$", candidate):
                # avoid operon-like acrAB
                if strain and candidate != strain:
                    substrain = candidate
                    query = (query[:m.start()] + query[m.end():]).strip()

    return strain, substrain, query


def _parse_deletions(query: str) -> List[Union[str, Tuple[str, str]]]:
    """Parse deletions from:
      - Δgene, delta gene, del-gene
      - ranges: Δ(fimB-fimE) -> ('fimB','fimE')
      - KO insertions: tolC::FRT (treat as tolC deletion)
      - operon shorthand with trailing '-' : acrAB- -> ['acrA','acrB']

    returns
    =======
    list
        Gene names.

    """
    deletions: List[Union[str, Tuple[str, str]]] = []

    # KO insertions "gene::something"
    for m in re.finditer(
        r"\b([A-Za-z][A-Za-z0-9._-]{1,})::[A-Za-z0-9._-]{2,}\b", 
        query,
    ):
        deletions.append(m.group(1))

    # Explicit Δ / delta / del
    for m in re.finditer(
        r"(?:Δ|delta|del)[\s-]*\(?([A-Za-z0-9_]+)(?:-([A-Za-z0-9._-]+))?\)?",
        query, 
        flags=re.IGNORECASE,
    ):
        a = _strip_punctuation(m.group(1))
        b = _strip_punctuation(m.group(2)) if m.group(2) else None
        if b:
            deletions.append((a, b))
        else:
            # Handle operon shorthand inside Δ token (e.g., ΔacrAB-)
            genes = _split_operon(a)
            deletions.extend(genes)

    # Standalone operon shorthand with trailing '-' (e.g., 'acrAB-')
    for m in re.finditer(
        r"\b([a-z]+[A-Z]{2,})-?\b",
        query,
    ):
        genes = _split_operon(m.group(1))
        deletions.extend(genes)

    normed = []
    for d in set(deletions):
        if isinstance(d, tuple):
            a, b = d
            normed.append((_normalize_gene(a), _normalize_gene(b)))
        else:
            normed.append(_normalize_gene(d))

    return sorted(set(normed), key=lambda x: x[0] if isinstance(x, tuple) else x)


def _normalize_gene(query: str) -> str:
    return query.strip().rstrip("-_,.;")


def _parse_mutations(
    query: str, 
    exclude: List[str]
) -> List[str]:
    """Very lightweight mutation token grabber, e.g. 'gyrA96', 'rpoB_S531L'.

    Excludes any tokens already classified as deletions.
    """
    mutations = []
    for m in re.finditer(
        r"\b([A-Za-z][A-Za-z0-9_]*\d+[A-Z0-9_]*)\b", 
        query,
    ):
        candidate = _strip_punctuation(m.group(1))
        if candidate not in exclude:
            mutations.append(candidate)
    return sorted(set(mutations))


# --- Public API ---

def parse_strain_label(
    query: Union[str, int]
) -> Strain:
    """Parse a free text strain name into components.

    Parameters
    ==========
    query : str
        Free text strain name.

    Returns
    =======
    Strain

    Examples
    ========
    >>> parse_strain_label("E. coli K-12 substr. MG1655 gyrA96 acrAB- Δ(fimB-fimE) ΔompF tolC::FRT")
    Strain(query='E. coli K-12 substr. MG1655 gyrA96 acrAB- Δ(fimB-fimE) ΔompF tolC::FRT', species='Escherichia coli', remainder='gyrA96 acrAB- Δ(fimB-fimE) ΔompF tolC::FRT', strain='K-12', substrain='MG1655', deletions=['acrA', 'acrB', ('fimB', 'fimE'), 'ompF', 'tolC'], mutations=['gyrA96'])
    >>> parse_strain_label("S. enterica serovar Typhimurium")
    Strain(query='S. enterica serovar Typhimurium', species='Salmonella enterica', remainder='', strain='Typhimurium', substrain=None, deletions=[], mutations=[])
    >>> parse_strain_label("S. aureus RN4220 Δspa")
    Strain(query='S. aureus RN4220 Δspa', species='Staphylococcus aureus', remainder='Δspa', strain='RN4220', substrain=None, deletions=['spa'], mutations=[])
    >>> parse_strain_label(83332)
    Strain(query='Mycobacterium tuberculosis H37Rv', species='Mycobacterium tuberculosis', remainder='', strain='H37Rv', substrain=None, deletions=[], mutations=[])

    """
    if isinstance(query, int) or (isinstance(query, str) and query.isdigit()):
        query = name_to_taxon_ncbi(query, key="sci_name")
    species, remainder = _extract_species(query, normalize=True)
    strain, substrain, remainder = _extract_strain_and_substrain(remainder)
    deletions = _parse_deletions(remainder)
    to_exclude = set()
    for d in deletions:
        if isinstance(d, tuple):
            to_exclude.update(d)
        else:
            to_exclude.add(d)
    mutations = _parse_mutations(
        remainder, 
        exclude=list(to_exclude),
    )

    return Strain(
        query=query,
        remainder=remainder,
        species=species,
        strain=strain,
        substrain=substrain,
        deletions=deletions,
        mutations=mutations,
    )
