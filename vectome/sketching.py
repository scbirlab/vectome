""""""

from typing import Optional
from functools import cache
import os

from carabiner import print_err
from sourmash import load_one_signature, MinHash, SourmashSignature, save_signatures

from .caching import CACHE_DIR
from .data import APPDATA_DIR
from .genomes import fetch_landmarks

DEFAULT_K: int = 101
DEFAULT_N: int = 200_000

@cache
def sketch_genome(
    file: str,
    k: int = DEFAULT_K,
    n: int = DEFAULT_N,
    force: bool = False,
    cache_dir: Optional[str] = None,
    _landmark: bool = False,  # prevents cache hits on landmark downloads
    **kwargs
):
    
    cache_dir = os.path.join(cache_dir or CACHE_DIR, "sketches")
    sketch_file = os.path.join(cache_dir, f"{os.path.basename(file)}_{n=}_{k=}.sig")

    if os.path.exists(sketch_file) and not force:
        print_err(f"Loading cached signature for {file} at {sketch_file}...", end=" ")
        try:
            mh = load_one_signature(sketch_file).minhash
        except ValueError: # no signatures to load
            print_err("failed!! Falling back to generating a sketch")
            return sketch_genome(
                file=file,
                k=k,
                n=n,
                force=True,
                cache_dir=os.path.dirname(cache_dir),
                _landmark=_landmark,
                **kwargs,
            )
        else:
            print_err("ok")
    else:
        from bioino import FastaCollection

        mh = MinHash(n=n, ksize=k, **kwargs)
        fasta = FastaCollection.from_file(file)
        for seq in fasta.sequences:
            mh.add_sequence(seq.sequence, force=True)

        sig = SourmashSignature(mh, name=os.path.basename(file))

        os.makedirs(os.path.dirname(sketch_file), exist_ok=True)
        print_err(f"Caching signature for {file} at {sketch_file}...", end=" ")
        with open(sketch_file, "w") as f:
            save_signatures([sig], f)
            print_err("ok")

    return mh


def sketch_landmarks(
    group: int = 0,
    check_spelling: bool = False,
    force: bool = False,
    cache_dir: Optional[str] = None
):
    from tqdm.auto import tqdm

    cache_dir = cache_dir or APPDATA_DIR
    landmark_info = fetch_landmarks(
        check_spelling=check_spelling,
        group=group,
        force=force,
        cache_dir=cache_dir,
    )

    landmarks = [
        info["files"]["fasta"] 
        for info in landmark_info
    ]
    return tuple(
        sketch_genome(
            file=f,
            force=force,
            cache_dir=cache_dir,
            _landmark=True,
        ) for f in tqdm(landmarks, desc="Sketching landmarks")
    )
