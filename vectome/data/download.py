"""Download pre-build landmark caches."""

from typing import Optional
from glob import glob
import tarfile
import os

from ..caching import CACHE_DIR
from ..http import download_url
from .. import __version__
from .manifest import (
    ASSET_NAME, 
    CHECKSUMS_NAME,
    DEFAULT_BASE_URL, 
    ENV_BASE_URL,
)

def _base_url() -> str:
    return os.environ.get(ENV_BASE_URL, DEFAULT_BASE_URL).rstrip("/")

def _release_url(
    version: str = __version__, 
    filename: str = ASSET_NAME
) -> str:
    return f"{_base_url()}/{version}/{filename}"

def download_landmark_cache(
    version: str = __version__,
    cache_dir: Optional[str] = None
) -> str:
    cache_dir = cache_dir or CACHE_DIR
    landmark_version = os.environ.get('VECTOME_LANDMARKS_VERSION', f"v{version}")
    dl_dir = os.path.join(cache_dir, "landmark-dl", version)
    dl_dir_temp = os.path.join(dl_dir, "temp")
    os.makedirs(dl_dir_temp, exist_ok=True)
    
    archive = download_url(
        _release_url(landmark_version),
        destination=dl_dir_temp,
    )

    with tarfile.open(archive, "r:*") as tf:
        tf.extractall(dl_dir)
    os.remove(archive)
    os.rmdir(dl_dir_temp)
    landmark_destination = os.path.join(cache_dir, "landmarks", landmark_version)
    os.makedirs(landmark_destination, exist_ok=True)
    os.rename(
        os.path.join(dl_dir, "landmarks"),
        landmark_destination,
    )
    sketch_destination = os.path.join(cache_dir, "sketches")
    os.makedirs(sketch_destination, exist_ok=True)
    for sketch_file in glob(os.path.join(dl_dir, "sketches", "*.sig")):
        os.rename(
            sketch_file, 
            os.path.join(sketch_destination, os.path.basename(sketch_file)),
        )
    os.rmdir(os.path.join(dl_dir, "sketches"))
    return cache_dir
    
