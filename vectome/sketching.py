""""""

from typing import Tuple
from functools import cache, partial
from io import StringIO
import os
from tempfile import NamedTemporaryFile

from carabiner import print_err
from sourmash import load_one_signature, MinHash, SourmashSignature, save_signatures

from .caching import CACHE_DIR, cache_lock
from .genomes import fetch_landmarks

DEFAULT_K: int = 21
DEFAULT_N: int = 10_000


def _load_sketch(
    file: str,
    sketch_file: str,
    quiet: bool = False
):
    if not quiet:
        print_err(f"Loading cached signature for {file} at {sketch_file}...", end=" ")
    mh = load_one_signature(sketch_file).minhash
    if not quiet:
        print_err("ok")
    return mh 


def _generate_sketch(
    file: str,
    k: int,
    n: int,
    **kwargs
):
    import gzip
    from bioino import FastaCollection

    mh = MinHash(n=n, ksize=k, **kwargs)
    opener = partial(gzip.open, mode="rb") if file.endswith(".gz") else partial(open, mode="r")
    with opener(file) as f:
        try:
            contents = f.read()
        except gzip.BadGzipFile as e:
            print_err(f"File '{file}' is not GZIP or is corrupted.")
            raise e
    if isinstance(contents, bytes):
        contents = contents.decode()
    fasta = FastaCollection.from_file(StringIO(contents))
    for seq in fasta.sequences:
        mh.add_sequence(seq.sequence, force=True)
    return mh


@cache
def sketch_genome(
    file: str,
    k: int = DEFAULT_K,
    n: int = DEFAULT_N,
    force: bool = False,
    quiet: bool = False,
    cache_dir: str = CACHE_DIR,
    _landmark: bool = False,  # prevents cache hits on landmark downloads
    **kwargs
) -> MinHash:
    from joblib import hash as joblib_hash
    sketch_id = joblib_hash((
        os.path.basename(file),
        n,
        k,
        kwargs,
    ))[:12]
    
    sketch_dir = os.path.join(cache_dir, "sketches")
    sketch_file = os.path.join(sketch_dir, f"{os.path.basename(file)}_{n=}_{k=}_{sketch_id}.sig")

    if os.path.exists(sketch_file) and not force:
        try:
            return _load_sketch(file, sketch_file, quiet=quiet)
        except ValueError:
            if not quiet:
                print_err(f"[WARN] Cached signature for {file} at {sketch_file} invalid; regenerating.")
    
    os.makedirs(sketch_dir, exist_ok=True)

    with cache_lock(
        key=(
            "sketch",
            os.path.basename(file),
            n,
            k,
            kwargs,
        ),
        cache_dir=cache_dir,
    ):
        # The process ahead of us may have populated/repaired it.
        if os.path.exists(sketch_file) and not force:
            try:
                return _load_sketch(file, sketch_file, quiet=quiet)
            except ValueError:
                pass

        mh = _generate_sketch(
            file,
            n=n,
            k=k,
            **kwargs,
        )
        sig = SourmashSignature(mh, name=os.path.basename(file))

        if not quiet:
            print_err(f"[INFO] Caching signature for {file} at {sketch_file}...", end=" ")
        
        with NamedTemporaryFile(
            mode="w",
            dir=sketch_dir,
            prefix=".tmp-",
            suffix=".sig",
            delete=False,
        ) as f:
            save_signatures([sig], f)
            os.replace(
                f.name,
                sketch_file,
            )
        if os.path.exists(f.name):
            os.remove(f.name)
        if not quiet:
            print_err("ok")

    return mh


def sketch_landmarks(
    group: int = 0,
    check_spelling: bool = False,
    force: bool = False,
    validate_fasta: bool = True,
    cache_dir: str = CACHE_DIR,
    max_workers: int = 1,
    **kwargs
) -> Tuple[MinHash]:
    from tqdm.auto import tqdm
    from tqdm.contrib.concurrent import process_map

    landmark_info = fetch_landmarks(
        check_spelling=check_spelling,
        group=group,
        force=force,
        allow_missing_files=not validate_fasta,
        cache_dir=cache_dir,
    )

    landmarks = [
        info["files"]["fasta"] 
        for info in landmark_info
    ]
    fn = partial(
        sketch_genome, 
        force=force,
        cache_dir=cache_dir,
        _landmark=True,
        **kwargs,
    )
    return process_map(
        fn, 
        landmarks, 
        max_workers=max_workers, 
        desc="Sketching landmarks",
    )
