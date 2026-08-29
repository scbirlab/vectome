"""Fetching remote data."""

from typing import Any, Callable, Mapping, Optional
from functools import cache, wraps
from importlib.metadata import metadata
import os
import random
from tempfile import NamedTemporaryFile
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
    cache_dir: str = CACHE_DIR,
    cache_subdir: str | None = None,
    skip_cache: bool = False
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
        time.sleep(random.uniform(wait * .5, wait * 1.5))
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
        except (
            requests.exceptions.ConnectionError,
            requests.exceptions.Timeout,
        ) as error:
            retry = True
            r = None
        else:
            retry = r.status_code in {
                429,
                500,
                502,
                503,
                504,
            }
            error = None
            if not quiet:
                print_err(f"[INFO] ({_try + 1} / {max_tries}) Tried {r.url}... {r.status_code}", end="")
            if not retry:
                r.raise_for_status()
                if not quiet:
                    print_err(f"... ok")
                return f(query, r, *args, **kwargs)
            else:
                if not quiet:
                    print_err(f"... retrying")

        next_try = _try + 1
        
        if next_try < max_tries:
            print_err("")
            if r is not None:
                retry_after = r.headers.get("Retry-After")
                try:
                    next_wait = float(retry_after)
                except ValueError:
                    next_wait = wait * 2.
            else:
                retry_after = None
                next_wait = random.uniform(
                    wait,
                    wait * 2.,
                )
            if not quiet:
                print_err(f"[INFO] Tried {next_try} / {max_tries} times...", end=" ")
            time.sleep(next_wait)
            return api_call(
                query=query, 
                params=params,
                wait=min(wait * 2., 30),
                _try=next_try,
                *args, **kwargs
            )
        if error is not None:
            print_err("stopping!")
            raise error
        r.raise_for_status()

    endpoint = f"{f.__module__}.{f.__qualname__}.{url}"
    api_call.__module__ = f.__module__
    api_call.__name__ = f"{f.__name__}_api_call"
    api_call.__qualname__ = f"{f.__qualname__}_api_call"

    def cached_call(runtime_cache_dir):
        if cache_subdir is not None:
            runtime_cache_dir = os.path.join(
                runtime_cache_dir,
                cache_subdir,
            )
        return locked_cache(
            api_call,
            cache_dir=os.path.join(
                runtime_cache_dir,
                "api_calls",
            ),
            namespace=endpoint,
            ignore=["quiet", "wait", "_try"],
        )
    

    @wraps(f)
    def wrapper(
        *args,
        cache_dir=cache_dir,
        **kwargs,
    ):
        if skip_cache:
            return api_call(*args, **kwargs)
        return cached_call(cache_dir)(*args, **kwargs)
        
    return wrapper


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
        with NamedTemporaryFile(dir=cache_dir, delete=False) as f:
            f.write(r.content)
            os.replace(f.name, destination)
        if os.path.exists(f.name):
            os.remove(f.name)
        return destination

    return _download_url(query=(url, destination), quiet=quiet, wait=wait)
