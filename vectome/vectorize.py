"""Convert species names into vectors."""

from typing import Iterable, Optional, Union
from functools import cache
import hashlib
import os

from bioino import FastaCollection
from carabiner import cast, print_err
from sourmash import MinHash
from tqdm.auto import tqdm

from .caching import CACHE_DIR
from .ncbi import fetch_landmarks, name_or_taxon_to_genome_info
from .sketching import sketch_genome


def _mix_u64(x: int) -> int:
    """64-bit mix function (xorshift* / splitmix-like) to derive secondary hashes
    deterministically from a base 64-bit integer. 
    
    This is purely to decorrelate index/sign computations across num_hash_fns.

    Examples
    ========
    >>> # _mix_u64 is deterministic on a 64-bit domain
    >>> hex(_mix_u64(0))
    '0x0'
    >>> hex(_mix_u64(1))
    '0xb456bcfc34c2cb2c'
    >>> hex(_mix_u64(0x123456789ABCDEF0))
    '0x18b8c062f6f42398'

    """
    x &= 0xFFFFFFFFFFFFFFFF
    x ^= (x >> 33)
    x = (x * 0xff51afd7ed558ccd) & 0xFFFFFFFFFFFFFFFF
    x ^= (x >> 33)
    x = (x * 0xc4ceb9fe1a85ec53) & 0xFFFFFFFFFFFFFFFF
    x ^= (x >> 33)
    return x & 0xFFFFFFFFFFFFFFFF


def _bucket_index(h: int, dim: int, salt: int) -> int:
    """Deterministic bucket index in [0, dim).
    
    We combine the 64-bit hash 'h' with a small integer salt 'salt'
    (0..num_hash_fns-1), then reduce modulo 'dim'.

    Examples
    ========
    >>> # _bucket_index/_bucket_sign are deterministic and stable
    >>> h = _mix_u64(0x1234)
    >>> h, _bucket_index(h, 1024, 0), _bucket_sign(h, 0)
    (11969492833970939502, 635, 1)
    >>> _bucket_index(h, 1024, 1), _bucket_sign(h, 1)
    (580, -1)
    >>> _bucket_index(h, 10, 2), _bucket_sign(h, 2)
    (3, -1)

    """
    y = (h ^ ((salt + 1) * 0x9E3779B97F4A7C15)) & 0xFFFFFFFFFFFFFFFF
    return int(y % dim)


def _bucket_sign(h: int, salt: int) -> int:
    """Deterministic sign bit (+1/-1) derived from SHA1(h || salt).
    
    Using cryptographic hashing here makes the sign unbiased and stable
    across Python versions/platforms.
    """
    b = h.to_bytes(8, "little", signed=False) + salt.to_bytes(4, "little", signed=False)
    return 1 if (int(hashlib.sha1(b).hexdigest(), 16) & 1) else -1    


def _vectorize_landmark(
    query_mh: MinHash,
    group: int = 0,
    **kwargs
):
    landmarks = [
        info["files"]["fasta"] for info in fetch_landmarks(group=group)
    ]
    landmark_mh = [
        sketch_genome(
            file=f,
        ) for f in landmarks
    ]

    return [query_mh.similarity(lm) for lm in landmark_mh]


def _vectorize_countsketch(
    query_mh: MinHash,
    dim: int = None, 
    num_hash_fns: int = 3,
    **kwargs
):
    """Convert a sourmash MinHash (which holds a set of 64-bit hashes) into a fixed-length
    real-valued vector using CountSketch/feature hashing.

    For each 64-bit hash h in the MinHash set:
        for k in 0..num_hash_fns-1:
            h_k = mix(h ^ k)                  # derive independent-ish stream
            j   = bucket_index(h_k, dim, k)  # bucket 0..dim-1
            s   = bucket_sign(h_k, k)        # +1 or -1
            v[j] += s

    After processing all hashes, L2-normalize v.

    Notes:
    - Using num_hash_fns in {2,3,4} reduces variance vs single-bucket placement.
    - Larger dim reduces collisions; 1024–4096 is a good starting range.

    Examples
    ========
    >>> # _vectorize_countsketch on a tiny "MinHash-like" dummy:
    >>> class _DummyMH:  # exposes .hashes like sourmash.MinHash
    ...     def __init__(self, ints): self.hashes = {int(x): 1 for x in ints}
    >>> v = _vectorize_countsketch(_DummyMH([0x1234, 0xBEEF]), dim=16, num_hash_fns=3)
    >>> len(v), float((v @ v) ** 0.5)  # length and unit norm
    (16, 1.0)
    >>> # Expected nonzeros after normalization (see analysis):
    >>> round(v[4], 3), round(v[10], 3), round(v[13], 3), round(v[0], 3), round(v[11], 3)
    (-0.5, -0.5, 0.5, -0.5, 0.0)

    """
    import numpy as np

    dim = dim or 4096
    vector = np.zeros((dim,))

    for _hash in query_mh.hashes.keys():
        base = _hash & 0xFFFFFFFFFFFFFFFF
        for i in range(num_hash_fns):
            hk = _mix_u64(base ^ i)
            idx = _bucket_index(hk, dim, i)
            sign = _bucket_sign(hk, i)
            vector[idx] += float(sign)

    # L2 normalize to decouple vector length from sketch size
    vector_norm = np.linalg.norm(vector)
    if vector_norm > 0:
        vector /= vector_norm

    return vector


def vectorize(
    query: Union[str, int, Iterable[Union[str, int]]],
    k: int = 51,
    method: str = "countsketch",
    projection: Optional[int] = None,
    seed: int = 42,
    cache_dir: Optional[str] = None,
    **kwargs
):
    from joblib import Memory
    import numpy as np

    cache_dir = cache_dir or CACHE_DIR
    mem = Memory(location=os.path.join(cache_dir, "vectors"))

    genome_info = [
        name_or_taxon_to_genome_info(
            q, 
            cache_dir=cache_dir,
        ) for q in cast(query, to=list)
    ]

    query_mh = [
        sketch_genome(
            file=info["files"]["fasta"],
        ) for info in genome_info
    ]

    if method == "landmark":
        fn = mem.cache(_vectorize_landmark)
    elif method == "countsketch":
        fn = mem.cache(_vectorize_countsketch)
    else:
        raise ValueError(f"Vectorization {method=} is not implemented.")

    vectors = np.stack([
        fn(mh, **kwargs) for mh in tqdm(query_mh, desc="Vectorizing genomes")
    ], axis=0)

    if projection is None:
        return vectors
    else:
        from numpy.random import default_rng

        generator = default_rng(seed=seed)
        projector = generator.normal(
            scale=1. / np.sqrt(projection), 
            size=(vectors.shape[1], projection),
        )
        print_err("Projecting...")
        return vectors @ projector
