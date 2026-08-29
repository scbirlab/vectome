"""Fetching remote data."""

from typing import Any, Callable, Mapping, Optional
from functools import cache, wraps
from importlib.metadata import metadata
import os
import time

from carabiner import print_err, pprint_dict
from carabiner.decorators import decorator_with_params
import requests
from requests import Response

from . import app_name, __version__, __author__
from .caching import CACHE_DIR, locked_cache


@decorator_with_params
def api_get(
    f: Callable[[str, Response], ...],
    url: str,
    max_tries: int = 3,
    query_key: str | None = None,
    default_params: Mapping[str, ...] | None = None,
    cache_dir: str = CACHE_DIR
) -> Callable[[str | None, dict | None], None]:
    default_params = default_params or {}
    url0 = url

    headers = {
        "User-Agent": f"{app_name}/{__version__}",
        "From": (
            metadata(app_name)["Author-email"]
            .replace(__author__, "")
            .strip()
        ),
    }

    def api_call(
        query=None, 
        params=None,
        quiet: bool = False,
        wait: float = .2,
        _try: int = 0,
        *args, **kwargs
    ):
        time.sleep(wait)
        params = default_params | (params or {})
        url = url0
        if query_key is not None and query is not None:
            params = params | {query_key: query}
        elif query is not None:
            url = url.format(query=query)
        if not quiet:
            pprint_dict(
                params,
                message=f"Downloading from '{url}' with the following parameters"
            )
        try:
            r = requests.get(url, params=params, headers=headers)
        except requests.exceptions.ConnectionError as e:
            print_err(e)
            next_try = _try + 1
            if not quiet:
                print_err(f"[INFO] Tried {next_try} / {max_tries} times...", end=" ")
            if next_try < max_tries:
                print_err("")
                return api_call(
                    query=query, 
                    params=params,
                    wait=wait * 2,
                    _try=next_try,
                    *args, **kwargs
                )
            else:
                print_err("stopping!")
                raise e
        else:
            if not quiet:
                print_err(f"Trying {r.url}... {r.status_code}", end="")
            r.raise_for_status()
            if not quiet:
                print_err(f"... ok")
            return f(query, r, *args, **kwargs)

    if cache_dir is not None and isinstance(cache_dir, str):
        api_call = locked_cache(
            api_call, 
            cache_dir=os.path.join(cache_dir, "api_calls"),
        )
        
    return wraps(f)(cache(api_call))


def download_url(
    url: str,
    destination: str,
    max_tries: int = 3,
    wait: float = .2,
    quiet: bool = False,
    cache_dir: str = CACHE_DIR
) -> str:

    @api_get(
        url=url,
        max_tries=max_tries,
        cache_dir=cache_dir,
    )
    def _download_url(
        query,
        r: Response
    ):
        with open(destination, 'wb') as f:
            f.write(r.content)
        return destination

    return _download_url(query=(url, destination), quiet=quiet, wait=wait)
