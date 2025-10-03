"""Load package data."""

from typing import Any, Dict, Optional
import os

APPDATA_DIR = os.path.dirname(__file__)

def load_landmarks(
    cache_dir: Optional[str] = None
) -> Dict[str, Any]:

    import yaml

    cache_dir = cache_dir or APPDATA_DIR

    landmarks_path = os.path.join(
        APPDATA_DIR, 
        "landmarks.yml",
    )

    with open(landmarks_path, "r") as f:
         data = yaml.safe_load(f)

    info = {}
    for item in data["groups"]:
        to_append = []
        for q in item["queries"]:
            if q.startswith("file://"):
                pathname = os.path.join(APPDATA_DIR, q.split("file://")[-1])
                with open(pathname, "r") as f:
                    for line in f:
                        to_append.append(line.rstrip())
            else:
                to_append.append(q.rstrip())
        info[item["name"]] = to_append
    return info


def landmark_info(
    cache_dir: Optional[str] = None
):
    info = {
        key: len(value) if isinstance(value, (list, tuple)) else value
        for key, value in load_landmarks(cache_dir=cache_dir).items()
    }
    cache_dir = cache_dir or APPDATA_DIR
    cache_dir = os.path.join(cache_dir, "landmarks")

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
