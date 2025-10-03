"""Fetching remote data."""

from typing import Any, Callable, Mapping, Optional
from functools import cache, wraps
import time

from carabiner import print_err, pprint_dict
from carabiner.decorators import decorator_with_params
import requests
from requests import Response


@decorator_with_params
def api_get(
    f: Callable[[Response], Any],
    url: str,
    query_key: Optional[str] = None,
    default_params: Optional[Mapping[str, Any]] = None,
    cache_dir: Optional[str] = None
):
    default_params = default_params or {}
    url0 = url

    def api_call(query=None, params=None, *args, **kwargs):
        time.sleep(.2)
        params = default_params | (params or {})
        url = url0
        if query_key is not None and query is not None:
            params = params | {query_key: query}
        elif query is not None:
            url = url.format(query=query)
        pprint_dict(
            params,
            message=f"Downloading from {url} with the following parameters"
        )
        r = requests.get(url, params=params)
        print_err(f"Trying {r.url}", end="")
        r.raise_for_status()
        print_err(f"... {r.status_code} ok")
        return f(query, r, *args, **kwargs)

    if cache_dir is None or not isinstance(cache_dir, str):
        return wraps(f)(cache(api_call))
    else:
        from joblib import Memory
        mem = Memory(location=cache_dir, verbose=0)
        return wraps(f)(cache(mem.cache(api_call)))
