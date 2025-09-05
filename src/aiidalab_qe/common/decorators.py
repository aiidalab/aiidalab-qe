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

    Usage
    -----
    Methods / free functions:
    >>> @cache_per_thread
    >>> def f(...): ...

    Properties (must use this order):
    >>> @cache_per_thread(invalidator="uuid")
    >>> @property
    >>> def f(...): ...
    """

    def cached_call(func, *args, **kwargs):
        # Get or initialize the cache for this thread
        cache = getattr(_thread_local, "cache", None)
        if cache is None:
            cache = {}
            _thread_local.cache = cache

        # Begin constructing the cache key
        key_parts = [func]

        if args:
            # Check if func is a method of a class
            if hasattr(args[0], "__dict__"):
                self = args[0]
                key_parts.append(id(self))

                # Include the invalidator attribute if specified
                if invalidator is not None:
                    key_parts.append(getattr(self, invalidator))

        # Include the arguments in the cache key
        if args or kwargs:
            key_parts.append(args)
            key_parts.append(frozenset(kwargs.items()))

        key = tuple(key_parts)

        # Check if the result is already cached
        if key not in cache:
            cache[key] = func(*args, **kwargs)

        return cache[key]

    def decorator(func):
        # Property case
        if isinstance(func, property):
            fget = func.fget
            if fget is None:
                raise TypeError("Property has no getter to wrap")

            @functools.wraps(fget)
            def wrapper(self):  # type: ignore
                return cached_call(fget, self)

            return property(wrapper)

        # method/function case
        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            return cached_call(func, *args, **kwargs)

        return wrapper

    return decorator
