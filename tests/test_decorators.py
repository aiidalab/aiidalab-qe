from threading import Thread

from aiidalab_qe.common.decorators import cache_per_thread


class DummyClass:
    def __init__(self, value):
        self.value = value

    @cache_per_thread("value")
    def compute(self, x):
        return self.value + x


def test_cache_per_thread_single_thread():
    obj = DummyClass(10)
    assert obj.compute(5) == 15  # First call, computes the result
    assert obj.compute(5) == 15  # Cached result, no recomputation


def test_cache_per_thread_different_args():
    obj = DummyClass(10)
    assert obj.compute(5) == 15  # First call with x=5
    assert obj.compute(3) == 13  # First call with x=3
    assert obj.compute(5) == 15  # Cached result for x=5


def test_cache_per_thread_different_threads():
    obj = DummyClass(300)  # must be > 256 to avoid Python integer interning
    results = []

    def thread_func(x):
        results.append(obj.compute(x))

    thread1 = Thread(target=thread_func, args=(5,))
    thread2 = Thread(target=thread_func, args=(5,))

    thread1.start()
    thread2.start()
    thread1.join()
    thread2.join()

    # Each thread should compute independently
    assert results[0] == 305
    assert results[1] == 305

    # Ensure results are not shared between threads
    assert results[0] is not results[1]


def test_cache_invalidation():
    obj = DummyClass(10)
    assert obj.compute(5) == 15  # First call, computes the result
    assert obj.compute(5) == 15  # Cached result, no recomputation
    obj.value = 20  # Invalidate cache by changing the invalidator attribute
    assert obj.compute(5) == 25  # Recomputes the result
