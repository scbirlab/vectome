"""Test single-flight caching."""

import multiprocessing as mp
import os
import time

import pytest

from vectome.caching import locked_cache


def _increment_counter(counter_file):
    """Increment a file-backed counter once.

    Deliberately slow enough that concurrent uncached calls overlap.
    """
    time.sleep(.1)

    try:
        with open(counter_file) as f:
            value = int(f.read())
    except FileNotFoundError:
        value = 0

    with open(counter_file, "w") as f:
        f.write(str(value + 1))

    return value + 1


def _call_cached(cache_dir, counter_file, queue):
    cached = locked_cache(
        _increment_counter,
        cache_dir=cache_dir,
    )
    queue.put(cached(counter_file))


def test_locked_cache_reuses_result(tmp_path):
    counter_file = tmp_path / "counter.txt"
    cache_dir = tmp_path / "cache"

    cached = locked_cache(
        _increment_counter,
        cache_dir=str(cache_dir),
    )

    first = cached(str(counter_file))
    second = cached(str(counter_file))

    assert first == 1
    assert second == 1
    assert counter_file.read_text() == "1"


def test_locked_cache_different_arguments_are_independent(tmp_path):
    cache_dir = tmp_path / "cache"

    counter_a = tmp_path / "a.txt"
    counter_b = tmp_path / "b.txt"

    cached = locked_cache(
        _increment_counter,
        cache_dir=str(cache_dir),
    )

    assert cached(str(counter_a)) == 1
    assert cached(str(counter_b)) == 1

    assert counter_a.read_text() == "1"
    assert counter_b.read_text() == "1"


def test_locked_cache_single_flight_across_processes(tmp_path):
    """
    Many simultaneous processes requesting the same uncached value
    should cause exactly one underlying function execution.
    """
    cache_dir = str(tmp_path / "cache")
    counter_file = str(tmp_path / "counter.txt")

    ctx = mp.get_context("spawn")
    queue = ctx.Queue()

    processes = [
        ctx.Process(
            target=_call_cached,
            args=(cache_dir, counter_file, queue),
        )
        for _ in range(32)
    ]

    for process in processes:
        process.start()

    for process in processes:
        process.join(timeout=30)
        assert process.exitcode == 0

    results = [queue.get(timeout=1) for _ in processes]

    assert results == [1] * len(processes)

    with open(counter_file) as f:
        assert f.read() == "1"


def test_locked_cache_allows_different_keys_concurrently(tmp_path):
    """
    Per-key locking must not become one global cache lock.
    """
    cache_dir = str(tmp_path / "cache")

    ctx = mp.get_context("spawn")
    queue = ctx.Queue()

    counter_files = [
        str(tmp_path / f"counter-{i}.txt")
        for i in range(4)
    ]

    start = time.monotonic()

    processes = [
        ctx.Process(
            target=_call_cached,
            args=(cache_dir, counter_file, queue),
        )
        for counter_file in counter_files
    ]

    for process in processes:
        process.start()

    for process in processes:
        process.join(timeout=30)
        assert process.exitcode == 0

    elapsed = time.monotonic() - start

    results = [queue.get(timeout=1) for _ in processes]

    assert results == [1] * len(processes)

    for counter_file in counter_files:
        with open(counter_file) as f:
            assert f.read() == "1"

    # Each underlying call sleeps for 0.1 s.
    # A single global lock would make four calls take ~0.4 s plus
    # process overhead. This is intentionally generous for CI.
    assert elapsed < 2.


def test_failed_computation_does_not_poison_lock(tmp_path):
    """
    flock should be released automatically when the wrapped function raises.
    """
    cache_dir = tmp_path / "cache"
    marker = tmp_path / "marker"

    def flaky(path):
        if not os.path.exists(path):
            open(path, "w").close()
            raise RuntimeError("first call fails")
        return 42

    cached = locked_cache(
        flaky,
        cache_dir=str(cache_dir),
    )

    with pytest.raises(RuntimeError, match="first call fails"):
        cached(str(marker))

    assert cached(str(marker)) == 42
