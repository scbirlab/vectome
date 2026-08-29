"""Load package data."""

from typing import Any, Dict, Optional
import os

from .. import __version__
from ..caching import CACHE_DIR

APPDATA_DIR = os.path.dirname(__file__)

def load_landmarks(
    cache_dir: str = APPDATA_DIR
) -> Dict[str, Any]:

    import yaml

    landmarks_path = os.path.join(
        cache_dir, 
        "landmarks.yml",
    )
    with open(landmarks_path, "r") as f:
         data = yaml.safe_load(f)

    info = {}
    for item in data["groups"]:
        to_append = []
        for q in item["queries"]:
            if q.startswith("file://"):
                pathname = os.path.join(cache_dir, q.split("file://")[-1])
                with open(pathname, "r") as f:
                    for line in f:
                        to_append.append(line.rstrip())
            else:
                to_append.append(q.rstrip())
        info[item["name"]] = to_append
    return info


def landmark_info(
    cache_dir: str = CACHE_DIR
):
    info = {
        key: len(value) if isinstance(value, (list, tuple)) else value
        for key, value in load_landmarks(cache_dir=APPDATA_DIR).items()
    }
    cache_dir = os.path.join(cache_dir, "landmarks", __version__)

    for key in info:
        group_cache = os.path.join(cache_dir, key)
        group_manifest = os.path.join(group_cache, "manifest.json")
        info[key] = {
            "landmarks": info[key],
            "manifest file": group_manifest,
            "built": os.path.exists(group_manifest),
        }
    info["meta"] = {
        "cache location": cache_dir,
        "cache exists": os.path.exists(cache_dir),
    }
    return info
