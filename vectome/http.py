"""Fetching remote data."""

from typing import Any, Callable, Mapping, Optional
from functools import cache, wraps
import os
import time

from carabiner import print_err, pprint_dict
from carabiner.decorators import decorator_with_params
import requests
from requests import Response


@decorator_with_params
def api_get(
    f: Callable[[str, Response], Any],
    url: str,
    max_tries: int = 3,
    query_key: Optional[str] = None,
    default_params: Optional[Mapping[str, Any]] = None,
    cache_dir: Optional[str] = None
) -> Callable[[Optional[str], Optional[dict]], Any]:
    default_params = default_params or {}
    url0 = url

    def api_call(
        query=None, 
        params=None,
        _try: int = 0,
        *args, **kwargs
    ):
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
        try:
            r = requests.get(url, params=params)
        except requests.exceptions.ConnectionError as e:
            next_try = _try + 1
            if next_try < max_tries:
                return api_call(
                    query=query, 
                    params=params,
                    _try=next_try,
                    *args, **kwargs
                )
            else:
                raise e
        else:
            print_err(f"Trying {r.url}", end="")
            r.raise_for_status()
            print_err(f"... {r.status_code} ok")
            return f(query, r, *args, **kwargs)

    if cache_dir is not None and isinstance(cache_dir, str):
        from joblib import Memory
        mem = Memory(location=os.path.join(cache_dir, "api_calls"), verbose=0)
        api_call = mem.cache(api_call)
        
    return wraps(f)(cache(api_call))
