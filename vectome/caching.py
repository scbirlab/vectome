"""Utilities for consistent caching."""

from collections.abc import Callable
from contextlib import contextmanager
import fcntl
import hashlib
from inspect import getsource
import os
from functools import wraps

from platformdirs import user_cache_dir

from . import app_name, __author__, __version__

CACHE_DIR = user_cache_dir(
    app_name, 
    __version__,
)


@contextmanager
def cache_lock(
    key: str,
    cache_dir: str = CACHE_DIR
):
    from joblib import hash as joblib_hash
    lock_dir = os.path.join(cache_dir, ".locks")
    os.makedirs(lock_dir, exist_ok=True)

    lock_id = joblib_hash(key)
    lock_file = os.path.join(lock_dir, f"{lock_id}.lock")

    with open(lock_file, "w") as f:
        fcntl.flock(f, fcntl.LOCK_EX)
        try:
            yield
        finally:
            fcntl.flock(f, fcntl.LOCK_UN)


def locked_cache(
    fn,
    cache_dir: str = CACHE_DIR,
    namespace: str | None = None,
    source: str | None = None,
    use_source: bool = False,
    **kwargs
):
    """Cache a function with single-flight locking across processes. 
    
    Concurrent calls with the same cache key are serialized while the first result is computed. 
    Calls with different keys can proceed independently. 
    
    Examples 
    ======== 
    >>> import tempfile 
    >>> calls = [] 
    >>> with tempfile.TemporaryDirectory() as cache_dir: 
    ...     def square(x): 
    ...         calls.append(x) 
    ...         return x * x 
    ... 
    ...     cached_square = locked_cache(square, cache_dir=cache_dir) 
    ...     cached_square(3), cached_square(3) 
    (9, 9)
    >>> len(calls)
    1
    
    Different arguments create different cache entries: 

    >>> with tempfile.TemporaryDirectory() as cache_dir: 
    ...     calls = [] 
    ...     def square(x): 
    ...         calls.append(x) 
    ...         return x * x 
    ... 
    ...     cached_square = locked_cache(square, cache_dir=cache_dir) 
    ...     cached_square(2), cached_square(4)
    (4, 16)
    >>> len(calls)
    2

    """
    from joblib import Memory

    memory = Memory(cache_dir, verbose=0, **kwargs)
    cached_fn = memory.cache(fn)

    namespace = namespace or (
        f"{fn.__module__}.{fn.__qualname__}"
    )
    if source is None:
        if all([
            use_source,
            not isinstance(fn, type(print)),
            isinstance(fn, Callable),
        ]):
            source = getsource(fn)
        else:
            source = ""
    
    @wraps(fn)
    def wrapper(*args, **kwargs):
        # Fast path: almost all calls after warmup.
        if cached_fn.check_call_in_cache(*args, **kwargs):
            return cached_fn(*args, **kwargs)

        with cache_lock(
            key=(namespace, source, args, kwargs),
            cache_dir=cache_dir, 
        ):
            # Crucial second check after acquiring the lock.
            if cached_fn.check_call_in_cache(*args, **kwargs):
                return cached_fn(*args, **kwargs)

            return cached_fn(*args, **kwargs)

    return wrapper
