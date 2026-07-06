#include <features/matching_gpu.h>

#include <stdexcept>

#ifdef OPENSFM_HAVE_CUDA

#include <features/matching_gpu_cuda.h>
#include <pybind11/pybind11.h>

#include <algorithm>
#include <climits>
#include <cstring>
#include <unordered_set>
#include <vector>

namespace features {

namespace {

using gpu::KNN2ResultInt;
using gpu::PairOffsets;

// Use at most 1/4 of device global memory for descriptor/result buffers.
constexpr size_t kMemoryReserveFraction = 4;

/// Apply Lowe's ratio test on integer Hamming distances.
/// For Hamming distances (linear, not squared), the test is:
///   best_dist < ratio * second_dist
std::vector<std::pair<int, int>> ApplyRatioTestHamming(const KNN2ResultInt& knn,
                                                       float lowes_ratio,
                                                       int n1,
                                                       int query_offset = 0,
                                                       int read_offset = 0) {
  std::vector<std::pair<int, int>> matches;
  matches.reserve(n1 / 4);
  for (int i = 0; i < n1; ++i) {
    int ri = read_offset + i;
    if (knn.best_idx[ri] >= 0 &&
        (float)knn.best_dist[ri] < lowes_ratio * (float)knn.second_dist[ri]) {
      matches.emplace_back(query_offset + i, knn.best_idx[ri]);
    }
  }
  return matches;
}

/// Convert a vector of (i,j) pairs to a [M x 2] numpy int array.
py::array_t<int> PairsToArray(const std::vector<std::pair<int, int>>& pairs) {
  py::array_t<int> result({(int)pairs.size(), 2});
  auto r = result.mutable_unchecked<2>();
  for (size_t k = 0; k < pairs.size(); ++k) {
    r(k, 0) = pairs[k].first;
    r(k, 1) = pairs[k].second;
  }
  return result;
}

/// Compute how many query descriptors we can process in one batch given
/// the device memory and the fixed f2 buffer size.
int ComputeHammingBatchSize(int device_idx, int n2, int n_words) {
  size_t avail = gpu::CudaDeviceGlobalMem(device_idx) / kMemoryReserveFraction;
  size_t f2_bytes = (size_t)n2 * n_words * sizeof(uint32_t);
  size_t per_query = (size_t)n_words * sizeof(uint32_t) + 3 * sizeof(int);
  if (avail <= f2_bytes) {
    return 1024;
  }
  size_t batch = (avail - f2_bytes) / per_query;
  return std::max(1024, (int)std::min(batch, (size_t)INT_MAX));
}

void CheckDeviceIndex(int device_idx) {
  if (!gpu::CudaAvailable()) {
    throw std::runtime_error("CUDA is not available");
  }
  if (device_idx < 0 || device_idx >= gpu::CudaNumDevices()) {
    throw std::runtime_error("Invalid CUDA device index");
  }
}

struct PairHash {
  size_t operator()(const std::pair<int, int>& p) const {
    return std::hash<long long>()(((long long)p.first << 32) | p.second);
  }
};

}  // namespace

bool gpu_matching_available() { return gpu::CudaAvailable(); }

int gpu_num_devices() { return gpu::CudaNumDevices(); }

py::array_t<int> match_hamming_gpu(py::array_t<uint32_t> f1,
                                   py::array_t<uint32_t> f2, float lowes_ratio,
                                   int device_idx) {
  CheckDeviceIndex(device_idx);

  const int n1 = static_cast<int>(f1.shape(0));
  const int n2 = static_cast<int>(f2.shape(0));
  const int n_words = static_cast<int>(f1.shape(1));

  if (n1 == 0 || n2 == 0) {
    return py::array_t<int>(std::vector<int>{0, 2});
  }
  if (n_words != static_cast<int>(f2.shape(1))) {
    throw std::runtime_error(
        "Binary descriptor word counts must match between f1 and f2");
  }

  const uint32_t* f1_data = f1.data();
  const uint32_t* f2_data = f2.data();

  std::vector<std::pair<int, int>> all_matches;
  {
    py::gil_scoped_release release;

    int batch_size = ComputeHammingBatchSize(device_idx, n2, n_words);
    for (int offset = 0; offset < n1; offset += batch_size) {
      int count = std::min(batch_size, n1 - offset);
      auto knn =
          gpu::HammingKNN2(device_idx, f1_data + (size_t)offset * n_words,
                           count, f2_data, n2, n_words);
      auto batch = ApplyRatioTestHamming(knn, lowes_ratio, count, offset);
      all_matches.insert(all_matches.end(), batch.begin(), batch.end());
    }
  }

  return PairsToArray(all_matches);
}

py::array_t<int> match_hamming_gpu_symmetric(py::array_t<uint32_t> f1,
                                             py::array_t<uint32_t> f2,
                                             float lowes_ratio,
                                             int device_idx) {
  CheckDeviceIndex(device_idx);

  const int n1 = static_cast<int>(f1.shape(0));
  const int n2 = static_cast<int>(f2.shape(0));
  const int n_words = static_cast<int>(f1.shape(1));

  if (n1 == 0 || n2 == 0) {
    return py::array_t<int>(std::vector<int>{0, 2});
  }
  if (n_words != static_cast<int>(f2.shape(1))) {
    throw std::runtime_error(
        "Binary descriptor word counts must match between f1 and f2");
  }

  const uint32_t* f1_data = f1.data();
  const uint32_t* f2_data = f2.data();

  std::vector<std::pair<int, int>> symmetric;
  {
    py::gil_scoped_release release;

    // Forward pass: for each descriptor in f1, find NN in f2.
    int batch_fwd = ComputeHammingBatchSize(device_idx, n2, n_words);
    std::vector<std::pair<int, int>> matches_ij;
    for (int offset = 0; offset < n1; offset += batch_fwd) {
      int count = std::min(batch_fwd, n1 - offset);
      auto knn =
          gpu::HammingKNN2(device_idx, f1_data + (size_t)offset * n_words,
                           count, f2_data, n2, n_words);
      auto batch = ApplyRatioTestHamming(knn, lowes_ratio, count, offset);
      matches_ij.insert(matches_ij.end(), batch.begin(), batch.end());
    }

    // Reverse pass: for each descriptor in f2, find NN in f1.
    int batch_rev = ComputeHammingBatchSize(device_idx, n1, n_words);
    std::unordered_set<std::pair<int, int>, PairHash> reverse_set;
    for (int offset = 0; offset < n2; offset += batch_rev) {
      int count = std::min(batch_rev, n2 - offset);
      auto knn =
          gpu::HammingKNN2(device_idx, f2_data + (size_t)offset * n_words,
                           count, f1_data, n1, n_words);
      auto batch = ApplyRatioTestHamming(knn, lowes_ratio, count, offset);
      for (auto& p : batch) {
        // p = (j_offset, i_idx) → store as (i_idx, j_offset)
        reverse_set.emplace(p.second, p.first);
      }
    }

    // Intersect forward and reverse matches.
    symmetric.reserve(std::min(matches_ij.size(), reverse_set.size()));
    for (auto& p : matches_ij) {
      if (reverse_set.count(p)) {
        symmetric.push_back(p);
      }
    }
  }

  return PairsToArray(symmetric);
}

py::list match_hamming_gpu_batch_symmetric(py::list f1_list, py::list f2_list,
                                           float lowes_ratio, int device_idx) {
  CheckDeviceIndex(device_idx);

  int n_pairs = static_cast<int>(py::len(f1_list));
  if (n_pairs != static_cast<int>(py::len(f2_list))) {
    throw std::runtime_error("f1_list and f2_list must have the same length");
  }
  if (n_pairs == 0) {
    return py::list();
  }

  // Extract per-pair metadata.  Keep py::array_t handles alive so data
  // pointers remain valid throughout (even after GIL release — the
  // reference keeps the numpy array alive).
  struct PairArrays {
    py::array_t<uint32_t> a1, a2;
    const uint32_t* p1;
    const uint32_t* p2;
    int n1, n2;
  };
  std::vector<PairArrays> pair_info(n_pairs);

  int n_words = 0;
  for (int i = 0; i < n_pairs; ++i) {
    pair_info[i].a1 = f1_list[i].cast<py::array_t<uint32_t>>();
    pair_info[i].a2 = f2_list[i].cast<py::array_t<uint32_t>>();
    pair_info[i].n1 = static_cast<int>(pair_info[i].a1.shape(0));
    pair_info[i].n2 = static_cast<int>(pair_info[i].a2.shape(0));
    pair_info[i].p1 = pair_info[i].a1.data();
    pair_info[i].p2 = pair_info[i].a2.data();
    if (n_words == 0 && pair_info[i].n1 > 0) {
      n_words = static_cast<int>(pair_info[i].a1.shape(1));
    }
  }
  if (n_words == 0) {
    py::list result;
    for (int i = 0; i < n_pairs; ++i) {
      result.append(py::array_t<int>(std::vector<int>{0, 2}));
    }
    return result;
  }
  // All non-empty arrays must agree on the word count.
  for (int i = 0; i < n_pairs; ++i) {
    if ((pair_info[i].n1 > 0 &&
         static_cast<int>(pair_info[i].a1.shape(1)) != n_words) ||
        (pair_info[i].n2 > 0 &&
         static_cast<int>(pair_info[i].a2.shape(1)) != n_words)) {
      throw std::runtime_error(
          "All binary descriptor arrays must have the same word count");
    }
  }

  // Memory budget: use 1/kMemoryReserveFraction of device global memory.
  // Per HammingKNN2Batched call the GPU buffers are:
  //   queries = total_queries * n_words * 4
  //   refs    = total_refs    * n_words * 4
  //   results = total_queries * 3 * 4   (best_idx, best_dist, second_dist)
  //   block_info ≈ negligible
  // We run forward (q=f1, r=f2) then reverse (q=f2, r=f1) sequentially,
  // so peak = (sum_n1 + sum_n2) * desc_bytes + max(sum_n1, sum_n2) * 12.
  size_t avail = gpu::CudaDeviceGlobalMem(device_idx) / kMemoryReserveFraction;
  size_t desc_bytes = (size_t)n_words * sizeof(uint32_t);
  constexpr size_t kResultBytesPerQuery = 3 * sizeof(int);

  // Greedily partition pairs into chunks that fit in GPU memory.
  struct Chunk {
    int start, end;  // [start, end) into pair_info
    size_t sum_n1, sum_n2;
  };
  std::vector<Chunk> chunks;
  {
    int chunk_start = 0;
    size_t sum_n1 = 0, sum_n2 = 0;
    for (int i = 0; i < n_pairs; ++i) {
      size_t new_n1 = sum_n1 + pair_info[i].n1;
      size_t new_n2 = sum_n2 + pair_info[i].n2;
      size_t peak = (new_n1 + new_n2) * desc_bytes +
                    std::max(new_n1, new_n2) * kResultBytesPerQuery;
      if (peak > avail && i > chunk_start) {
        chunks.push_back({chunk_start, i, sum_n1, sum_n2});
        chunk_start = i;
        sum_n1 = pair_info[i].n1;
        sum_n2 = pair_info[i].n2;
      } else {
        sum_n1 = new_n1;
        sum_n2 = new_n2;
      }
    }
    if (chunk_start < n_pairs) {
      chunks.push_back({chunk_start, n_pairs, sum_n1, sum_n2});
    }
  }

  // Per-pair symmetric match results (populated under GIL release).
  std::vector<std::vector<std::pair<int, int>>> all_sym(n_pairs);

  {
    py::gil_scoped_release release;

    for (size_t ci = 0; ci < chunks.size(); ++ci) {
      auto& chunk = chunks[ci];
      int count = chunk.end - chunk.start;
      int total_n1 = static_cast<int>(chunk.sum_n1);
      int total_n2 = static_cast<int>(chunk.sum_n2);

      // Build chunk-local flat buffers and offsets.
      std::vector<uint32_t> chunk_f1((size_t)total_n1 * n_words);
      std::vector<uint32_t> chunk_f2((size_t)total_n2 * n_words);
      std::vector<PairOffsets> fwd_offsets(count);

      int off1 = 0, off2 = 0;
      for (int j = 0; j < count; ++j) {
        int i = chunk.start + j;
        fwd_offsets[j].f1_offset = off1;
        fwd_offsets[j].n1 = pair_info[i].n1;
        fwd_offsets[j].f2_offset = off2;
        fwd_offsets[j].n2 = pair_info[i].n2;
        if (pair_info[i].n1 > 0) {
          std::memcpy(chunk_f1.data() + (size_t)off1 * n_words, pair_info[i].p1,
                      (size_t)pair_info[i].n1 * n_words * sizeof(uint32_t));
        }
        if (pair_info[i].n2 > 0) {
          std::memcpy(chunk_f2.data() + (size_t)off2 * n_words, pair_info[i].p2,
                      (size_t)pair_info[i].n2 * n_words * sizeof(uint32_t));
        }
        off1 += pair_info[i].n1;
        off2 += pair_info[i].n2;
      }

      // Reverse offsets: query=f2, ref=f1.
      std::vector<PairOffsets> rev_offsets(count);
      for (int j = 0; j < count; ++j) {
        rev_offsets[j].f1_offset = fwd_offsets[j].f2_offset;
        rev_offsets[j].n1 = fwd_offsets[j].n2;
        rev_offsets[j].f2_offset = fwd_offsets[j].f1_offset;
        rev_offsets[j].n2 = fwd_offsets[j].n1;
      }

      // Forward: queries=f1, refs=f2.
      auto fwd = gpu::HammingKNN2Batched(device_idx, chunk_f1.data(), total_n1,
                                         chunk_f2.data(), total_n2,
                                         fwd_offsets, n_words);

      // Reverse: queries=f2, refs=f1.
      auto rev = gpu::HammingKNN2Batched(device_idx, chunk_f2.data(), total_n2,
                                         chunk_f1.data(), total_n1,
                                         rev_offsets, n_words);

      // Per-pair: ratio test + symmetric intersection.
      for (int j = 0; j < count; ++j) {
        int i = chunk.start + j;
        if (fwd_offsets[j].n1 == 0 || fwd_offsets[j].n2 == 0) {
          continue;
        }

        auto fwd_m = ApplyRatioTestHamming(fwd, lowes_ratio, fwd_offsets[j].n1,
                                           0, fwd_offsets[j].f1_offset);
        auto rev_m = ApplyRatioTestHamming(rev, lowes_ratio, rev_offsets[j].n1,
                                           0, rev_offsets[j].f1_offset);

        std::unordered_set<std::pair<int, int>, PairHash> rev_set;
        for (auto& p : rev_m) {
          rev_set.emplace(p.second, p.first);
        }

        auto& sym = all_sym[i];
        sym.reserve(std::min(fwd_m.size(), rev_set.size()));
        for (auto& p : fwd_m) {
          if (rev_set.count(p)) {
            sym.push_back(p);
          }
        }
      }
    }
  }
  // GIL re-acquired.

  py::list result;
  for (int i = 0; i < n_pairs; ++i) {
    result.append(PairsToArray(all_sym[i]));
  }
  return result;
}

}  // namespace features

#else  // !OPENSFM_HAVE_CUDA

namespace features {

bool gpu_matching_available() { return false; }

int gpu_num_devices() { return 0; }

py::array_t<int> match_hamming_gpu(py::array_t<uint32_t> /*f1*/,
                                   py::array_t<uint32_t> /*f2*/,
                                   float /*lowes_ratio*/, int /*device_idx*/) {
  throw std::runtime_error(
      "GPU matching is not available (built without CUDA support)");
}

py::array_t<int> match_hamming_gpu_symmetric(py::array_t<uint32_t> /*f1*/,
                                             py::array_t<uint32_t> /*f2*/,
                                             float /*lowes_ratio*/,
                                             int /*device_idx*/) {
  throw std::runtime_error(
      "GPU matching is not available (built without CUDA support)");
}

py::list match_hamming_gpu_batch_symmetric(py::list /*f1_list*/,
                                           py::list /*f2_list*/,
                                           float /*lowes_ratio*/,
                                           int /*device_idx*/) {
  throw std::runtime_error(
      "GPU matching is not available (built without CUDA support)");
}

}  // namespace features

#endif  // OPENSFM_HAVE_CUDA
