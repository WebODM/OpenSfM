#pragma once

#include <foundation/python_types.h>

#include <cstdint>

namespace features {

/// Returns true if at least one CUDA device is available for matching.
bool gpu_matching_available();

/// Returns the number of available CUDA devices.
int gpu_num_devices();

/// Brute-force Hamming matching on GPU (CUDA) for binary descriptors
/// with Lowe's ratio test.
///
/// Finds the nearest neighbour in f2 for each descriptor in f1, keeping
/// only matches that pass the ratio test (best_dist < ratio * second_dist).
///
/// @param f1          [N1 x W] uint32 packed binary descriptors (W words).
/// @param f2          [N2 x W] uint32 packed binary descriptors.
/// @param lowes_ratio Lowe's ratio threshold.
/// @param device_idx  CUDA device index (0 = first GPU).
/// @return [M x 2] int32 array of (query_idx, ref_idx) matches.
py::array_t<int> match_hamming_gpu(py::array_t<uint32_t> f1,
                                   py::array_t<uint32_t> f2, float lowes_ratio,
                                   int device_idx = 0);

/// Symmetric Hamming matching on GPU (CUDA): matches in both directions
/// and keeps only mutually consistent pairs.
///
/// @param f1          [N1 x W] uint32 packed binary descriptors.
/// @param f2          [N2 x W] uint32 packed binary descriptors.
/// @param lowes_ratio Lowe's ratio threshold.
/// @param device_idx  CUDA device index.
/// @return [M x 2] int32 array of (idx_in_f1, idx_in_f2) matches.
py::array_t<int> match_hamming_gpu_symmetric(py::array_t<uint32_t> f1,
                                             py::array_t<uint32_t> f2,
                                             float lowes_ratio,
                                             int device_idx = 0);

/// Batched symmetric Hamming matching on GPU (CUDA).
///
/// Matches all pairs in a single kernel dispatch for maximum GPU occupancy.
///
/// @param f1_list     List of [Ni x W] uint32 packed binary descriptor arrays.
/// @param f2_list     List of [Mi x W] uint32 packed binary descriptor arrays.
/// @param lowes_ratio Lowe's ratio threshold.
/// @param device_idx  CUDA device index.
/// @return List of [Ki x 2] int32 arrays of (idx_in_f1, idx_in_f2) matches.
py::list match_hamming_gpu_batch_symmetric(py::list f1_list, py::list f2_list,
                                           float lowes_ratio,
                                           int device_idx = 0);

}  // namespace features
