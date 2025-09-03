from __future__ import annotations

import functools
import threading

_thread_local = threading.local()


def cache_per_thread(invalidator: str | None = None):
    """Decorator for per-thread caching of functions, methods, or properties.

    Parameters
    ----------
    `invalidator` : `str | None`
        The name of the attribute to watch for changes.
        If the attribute changes, the cache will be invalidated.
    """

    def decorator(func):
        @functools.wraps(func)
        def wrapper(self, *args, **kwargs):
            cache = getattr(_thread_local, "cache", None)
            if cache is None:
                cache = {}
                _thread_local.cache = cache

            key_parts = [id(self), func]
            if invalidator is not None:
                key_parts.append(getattr(self, invalidator))
            if args or kwargs:
                key_parts.append(args)
                key_parts.append(frozenset(kwargs.items()))
            key = tuple(key_parts)

            if key not in cache:
                cache[key] = func(self, *args, **kwargs)

            return cache[key]

        if isinstance(func, property):
            return property(wrapper)

        return wrapper

    return decorator
