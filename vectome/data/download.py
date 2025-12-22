"""Download pre-build landmark caches."""

from typing import Optional
from glob import glob
import json
import tarfile
import os
import shutil

from requests import HTTPError

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
    suffix: str,
    version: str = __version__,
    filename: str = ASSET_NAME
) -> str:
    return f"{_base_url()}/{version}/{filename}-{suffix}.tar.gz"


def download_landmark_cache(
    suffix: str,
    version: str = __version__,
    quiet: bool = True,
    cache_dir: Optional[str] = None
) -> str:
    cache_dir = cache_dir or CACHE_DIR
    landmark_version = os.environ.get('VECTOME_LANDMARKS_VERSION', f"v{version}")
    dl_dir = os.path.join(cache_dir, "landmark-dl", landmark_version)
    dl_dir_temp = os.path.join(dl_dir, "temp")
    os.makedirs(dl_dir_temp, exist_ok=True)
    
    url = _release_url(suffix=suffix, version=landmark_version)
    try:
        archive = download_url(
            url,
            quiet=quiet,
            destination=os.path.join(dl_dir_temp, os.path.basename(url)),
        )
    except HTTPError as e:
        os.rmdir(dl_dir_temp)
        os.rmdir(dl_dir)
        raise e

    with tarfile.open(archive, "r:*") as tf:
        tf.extractall(dl_dir)
    os.remove(archive)
    os.rmdir(dl_dir_temp)
    landmark_destination = os.path.join(cache_dir, "landmarks", landmark_version)
    for landmark_dir in glob(os.path.join(dl_dir, "landmarks", "*", "group-*")):
        this_destination = os.path.join(landmark_destination, os.path.basename(landmark_dir))
        if os.path.exists(this_destination):
            shutil.rmtree(this_destination)
        shutil.move(
            landmark_dir,
            this_destination,
        )
        with open(os.path.join(this_destination, "manifest.json"), "r") as f:
            d = json.load(f)
        for item in d:
            for key in item["files"]:
                item["files"][key] = os.path.join(this_destination, os.path.basename(item["files"][key]))
        with open(os.path.join(this_destination, "manifest.json"), "w") as f:
            json.dump(d, f, indent=4)
    sketch_destination = os.path.join(cache_dir, "sketches")
    os.makedirs(sketch_destination, exist_ok=True)
    for sketch_file in glob(os.path.join(dl_dir, "sketches", "*.sig")):
        os.rename(
            sketch_file, 
            os.path.join(sketch_destination, os.path.basename(sketch_file)),
        )
    os.rmdir(os.path.join(dl_dir, "sketches"))
    return cache_dir
    
