"""Descriptor matching implemented with NumPy."""

# from __future__ import annotations

# import numpy as np


# def match_descriptors(
#     descriptors0: np.ndarray,
#     descriptors1: np.ndarray,
#     min_similarity: float = 0.82,
# ) -> np.ndarray:
#     """Return mutual nearest descriptor matches as ``(index0, index1)`` pairs.

#     Descriptors are expected to be row-wise L2-normalized.  The default
#     similarity threshold matches the official sparse XFeat implementation.
#     """
#     descriptors0 = np.asarray(descriptors0, dtype=np.float32)
#     descriptors1 = np.asarray(descriptors1, dtype=np.float32)
#     if descriptors0.ndim != 2 or descriptors1.ndim != 2:
#         raise ValueError("Descriptors must be two-dimensional arrays.")
#     if descriptors0.shape[1] != descriptors1.shape[1]:
#         raise ValueError("Descriptor dimensions must agree.")
#     if not -1.0 <= min_similarity <= 1.0:
#         raise ValueError("min_similarity must be in [-1, 1].")
#     if not len(descriptors0) or not len(descriptors1):
#         return np.empty((0, 2), dtype=np.int64)

#     similarity = descriptors0 @ descriptors1.T
#     nearest1 = np.argmax(similarity, axis=1)
#     nearest0 = np.argmax(similarity, axis=0)
#     indices0 = np.arange(len(descriptors0))
#     keep = (nearest0[nearest1] == indices0) & (
#         similarity[indices0, nearest1] >= min_similarity
#     )
#     return np.column_stack((indices0[keep], nearest1[keep])).astype(
#         np.int64, copy=False
#     )
