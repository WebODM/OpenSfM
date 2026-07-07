# pyre-strict
"""Tests for CUDA-accelerated brute-force Hamming descriptor matching."""

import numpy as np
import pytest
from numpy.typing import NDArray
from typing import List, Set, Tuple

from opensfm import pyfeatures


WORDS: int = 4  # 128-bit binary descriptors


def _random_binary_descriptors(n: int, seed: int) -> NDArray:
    rng = np.random.RandomState(seed)
    return rng.randint(0, 2**32, size=(n, WORDS), dtype=np.uint32)


def _plant_matches(
    f1: NDArray, f2: NDArray, n: int, offset1: int, offset2: int
) -> None:
    """Copy descriptors from f1[offset1:offset1+n] into f2[offset2:offset2+n]
    with one bit flipped (Hamming distance 1)."""
    f2[offset2: offset2 + n] = f1[offset1: offset1 + n]
    f2[offset2: offset2 + n, 0] ^= 1


def _cpu_hamming_matches(
    f1: NDArray, f2: NDArray, lowes_ratio: float
) -> Set[Tuple[int, int]]:
    """Brute-force Hamming 2-NN with Lowe's ratio test, for reference."""
    n1, n2 = len(f1), len(f2)
    xor = f1[:, None, :] ^ f2[None, :, :]
    dists = np.unpackbits(
        xor.view(np.uint8).reshape(n1, n2, -1), axis=2
    ).sum(axis=2)

    matches = set()
    for i in range(n1):
        order = np.argsort(dists[i], kind="stable")
        best, second = order[0], order[1]
        if dists[i][best] < lowes_ratio * dists[i][second]:
            matches.add((i, int(best)))
    return matches


def _to_set(result: NDArray) -> Set[Tuple[int, int]]:
    return {(int(r[0]), int(r[1])) for r in result}


gpu_available: bool = pyfeatures.gpu_matching_available()
skip_no_gpu = pytest.mark.skipif(not gpu_available, reason="CUDA not available")


@skip_no_gpu
def test_basic_matching() -> None:
    """Planted near-duplicate descriptors should be matched."""
    f1 = _random_binary_descriptors(200, 42)
    f2 = _random_binary_descriptors(300, 99)
    _plant_matches(f1, f2, 10, 0, 50)

    result = pyfeatures.match_hamming_gpu(f1, f2, 0.8)
    assert result.shape[1] == 2
    assert result.shape[0] >= 10

    matches_set = _to_set(result)
    for i in range(10):
        assert (i, 50 + i) in matches_set


@skip_no_gpu
def test_symmetric_matching() -> None:
    """Symmetric matching should return only mutually-consistent pairs."""
    f1 = _random_binary_descriptors(150, 11)
    f2 = _random_binary_descriptors(200, 22)
    _plant_matches(f1, f2, 5, 10, 30)

    sym_result = pyfeatures.match_hamming_gpu_symmetric(f1, f2, 0.8)
    asym_result = pyfeatures.match_hamming_gpu(f1, f2, 0.8)

    assert sym_result.shape[0] <= asym_result.shape[0]
    assert sym_result.shape[0] >= 5

    sym_set = _to_set(sym_result)
    for i in range(5):
        assert (10 + i, 30 + i) in sym_set


@skip_no_gpu
def test_empty_input() -> None:
    """Empty descriptor arrays should return empty matches."""
    f1 = np.zeros((0, WORDS), dtype=np.uint32)
    f2 = _random_binary_descriptors(50, 1)

    result = pyfeatures.match_hamming_gpu(f1, f2, 0.8)
    assert result.shape[0] == 0
    assert result.shape[1] == 2

    result2 = pyfeatures.match_hamming_gpu_symmetric(f1, f2, 0.8)
    assert result2.shape[0] == 0


@skip_no_gpu
def test_consistency_with_cpu() -> None:
    """GPU matches should be identical to a CPU reference implementation."""
    f1 = _random_binary_descriptors(80, 44)
    f2 = _random_binary_descriptors(100, 55)
    _plant_matches(f1, f2, 8, 0, 10)

    gpu_set = _to_set(pyfeatures.match_hamming_gpu(f1, f2, 0.75))
    cpu_set = _cpu_hamming_matches(f1, f2, 0.75)

    assert gpu_set == cpu_set


@skip_no_gpu
def test_batch_symmetric_matches_single_pair() -> None:
    """Batched matching should equal per-pair symmetric matching."""
    sizes = [(120, 200), (0, 80), (350, 0), (64, 512)]
    f1_list: List[NDArray] = []
    f2_list: List[NDArray] = []
    for p, (n1, n2) in enumerate(sizes):
        f1 = _random_binary_descriptors(n1, 100 + p)
        f2 = _random_binary_descriptors(n2, 200 + p)
        if n1 >= 20 and n2 >= 40:
            _plant_matches(f1, f2, 10, 0, 20)
        f1_list.append(f1)
        f2_list.append(f2)

    batch_results = pyfeatures.match_hamming_gpu_batch_symmetric(
        f1_list, f2_list, 0.8
    )
    assert len(batch_results) == len(sizes)

    for p, (n1, n2) in enumerate(sizes):
        if n1 == 0 or n2 == 0:
            assert batch_results[p].shape[0] == 0
            continue
        single = pyfeatures.match_hamming_gpu_symmetric(
            f1_list[p], f2_list[p], 0.8
        )
        assert _to_set(batch_results[p]) == _to_set(single)
