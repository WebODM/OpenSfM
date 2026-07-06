#pragma once

// Internal interface between the host-side matching logic (matching_gpu.cc)
// and the CUDA implementation (matching_gpu_cuda.cu).  Deliberately free of
// CUDA and Python types so matching_gpu.cc compiles with a plain C++
// compiler and matching_gpu_cuda.cu never touches the Python C API.
//
// These functions are only linked in when OPENSFM_HAVE_CUDA is defined;
// matching_gpu.cc must not reference them otherwise.

#include <cstddef>
#include <cstdint>
#include <vector>

namespace features {
namespace gpu {

/// Per-pair descriptor offsets for the batched kernel.  Offsets are in
/// descriptors (rows), not words.
struct PairOffsets {
  int f1_offset;
  int n1;
  int f2_offset;
  int n2;
};

/// 2-nearest-neighbour result with integer (Hamming) distances.
struct KNN2ResultInt {
  std::vector<int> best_idx;
  std::vector<int> best_dist;
  std::vector<int> second_dist;
};

/// Returns true if at least one usable CUDA device was found (driver
/// present, device compute capability supported by the compiled kernels).
bool CudaAvailable();

/// Number of usable CUDA devices.
int CudaNumDevices();

/// Total global memory of the given usable device, in bytes.
size_t CudaDeviceGlobalMem(int device_idx);

/// For each descriptor in f1, find its 2 nearest neighbours in f2 by
/// Hamming distance.  Blocking; throws std::runtime_error on CUDA errors.
/// Must be called without the Python GIL held.
KNN2ResultInt HammingKNN2(int device_idx, const uint32_t* f1, int n1,
                          const uint32_t* f2, int n2, int n_words);

/// Batched KNN2 across many (f1, f2) pairs in a single kernel dispatch.
/// queries/refs are the concatenated descriptors of all pairs; `pairs`
/// describes each pair's slice.  Results are indexed by absolute query row.
/// Must be called without the Python GIL held.
KNN2ResultInt HammingKNN2Batched(int device_idx, const uint32_t* queries,
                                 int total_queries, const uint32_t* refs,
                                 int total_refs,
                                 const std::vector<PairOffsets>& pairs,
                                 int n_words);

}  // namespace gpu
}  // namespace features
