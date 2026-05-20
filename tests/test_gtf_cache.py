"""LRU bound on the GTF parse cache.

The cache must:
1. Re-use parsed GTFs when keyed by the same canonical path (hit).
2. Evict the least-recently-used entry once size exceeds the limit.
3. Refresh recency on hit (so a touched entry is not the next to evict).
4. Resize on demand, evicting eagerly if the new max is smaller.
"""

from __future__ import annotations

import pytest

from epykit import annotate as A


@pytest.fixture(autouse=True)
def _isolate_cache():
    """Snapshot and restore module-level cache state around each test."""
    saved_items = list(A._GTF_CACHE.items())
    saved_size = A._GTF_CACHE_MAX_SIZE
    A._GTF_CACHE.clear()
    try:
        yield
    finally:
        A._GTF_CACHE.clear()
        A._GTF_CACHE.update(saved_items)
        A.set_gtf_cache_size(saved_size)


def test_lru_evicts_oldest_when_full():
    A.set_gtf_cache_size(2)
    A._gtf_cache_put("a", ("genes_a", "exons_a"))
    A._gtf_cache_put("b", ("genes_b", "exons_b"))
    A._gtf_cache_put("c", ("genes_c", "exons_c"))
    assert "a" not in A._GTF_CACHE
    assert "b" in A._GTF_CACHE
    assert "c" in A._GTF_CACHE


def test_lru_hit_refreshes_recency():
    A.set_gtf_cache_size(2)
    A._gtf_cache_put("a", ("genes_a", "exons_a"))
    A._gtf_cache_put("b", ("genes_b", "exons_b"))
    # Touch 'a' so 'b' becomes the LRU entry.
    assert A._gtf_cache_get("a") == ("genes_a", "exons_a")
    A._gtf_cache_put("c", ("genes_c", "exons_c"))
    assert "b" not in A._GTF_CACHE
    assert "a" in A._GTF_CACHE
    assert "c" in A._GTF_CACHE


def test_resize_evicts_eagerly():
    A.set_gtf_cache_size(3)
    A._gtf_cache_put("a", ("genes_a", "exons_a"))
    A._gtf_cache_put("b", ("genes_b", "exons_b"))
    A._gtf_cache_put("c", ("genes_c", "exons_c"))
    A.set_gtf_cache_size(1)
    assert len(A._GTF_CACHE) == 1
    assert "c" in A._GTF_CACHE  # most-recently inserted survives


def test_set_size_rejects_zero():
    with pytest.raises(ValueError):
        A.set_gtf_cache_size(0)


def test_get_miss_returns_none():
    assert A._gtf_cache_get("nonexistent") is None


def test_put_existing_key_updates_value_and_recency():
    A.set_gtf_cache_size(2)
    A._gtf_cache_put("a", ("genes_a1", "exons_a1"))
    A._gtf_cache_put("b", ("genes_b", "exons_b"))
    A._gtf_cache_put("a", ("genes_a2", "exons_a2"))  # update
    assert A._gtf_cache_get("a") == ("genes_a2", "exons_a2")
    # 'b' should now be LRU; inserting 'c' evicts it.
    A._gtf_cache_put("c", ("genes_c", "exons_c"))
    assert "b" not in A._GTF_CACHE
    assert "a" in A._GTF_CACHE
