// CUDA implementation of brute-force Hamming descriptor matching.
//
// Compiled only when OPENSFM_HAVE_CUDA is defined (see features
// CMakeLists.txt).  Kernels are built for compute capability 5.0
// (NVIDIA Maxwell, e.g. GeForce GTX 950M) plus PTX so newer GPUs can
// JIT-compile.  No Python types here — the interface with the host
// logic is matching_gpu_cuda.h.

#include <cuda_runtime.h>
#include <features/matching_gpu_cuda.h>

#include <algorithm>
#include <iostream>
#include <mutex>
#include <stdexcept>
#include <string>
#include <vector>

namespace features {
namespace gpu {

namespace {

// Number of reference descriptors staged in shared memory per tile.
// 256 descriptors × 4 words × 4 bytes = 4 KB for 128-bit descriptors —
// small enough for any supported GPU (sm_50 has 48 KB per block).
constexpr int kHammingTileSize = 256;

// Threads per block; one thread handles one query descriptor.
constexpr int kBlockSize = 256;

// Maximum descriptor word count supported by the generic (non-templated)
// kernel: 32 words = 1024 bits.  The pipeline uses 4 words (128 bits).
constexpr int kMaxWords = 32;

// Oldest device architecture the embedded kernels can run on: the minimum
// of the architectures selected at build time (CMAKE_CUDA_ARCHITECTURES),
// as sm_XY * 10.  PTX for this architecture is embedded too, so any newer
// device can JIT-compile.
#ifdef __CUDA_ARCH_LIST__
constexpr int MinCompiledArch() {
  constexpr int archs[] = {__CUDA_ARCH_LIST__};
  int min_arch = archs[0];
  for (int arch : archs) {
    if (arch < min_arch) {
      min_arch = arch;
    }
  }
  return min_arch;
}
constexpr int kMinSupportedArch = MinCompiledArch();
#else
constexpr int kMinSupportedArch = 500;  // sm_50 (Maxwell, GTX 950M)
#endif

/// Check a CUDA API call result and throw on failure.
inline void CheckCuda(cudaError_t err, const char* context) {
  if (err != cudaSuccess) {
    throw std::runtime_error(std::string("CUDA error in ") + context + ": " +
                             cudaGetErrorString(err) + " (code " +
                             std::to_string((int)err) + ")");
  }
}

#define CUDA_CHECK(call) CheckCuda((call), #call)

/// RAII device buffer so error paths never leak GPU memory.
class DeviceBuffer {
 public:
  explicit DeviceBuffer(size_t bytes) {
    CheckCuda(cudaMalloc(&ptr_, bytes), "cudaMalloc");
  }
  ~DeviceBuffer() {
    if (ptr_ != nullptr) {
      cudaFree(ptr_);
    }
  }
  DeviceBuffer(const DeviceBuffer&) = delete;
  DeviceBuffer& operator=(const DeviceBuffer&) = delete;

  template <typename T>
  T* as() const {
    return static_cast<T*>(ptr_);
  }

 private:
  void* ptr_ = nullptr;
};

/// Discovers usable CUDA devices once per process.  A device is usable if
/// its compute capability is at least 5.0 — the oldest architecture the
/// kernels are compiled for.  Construction never throws: with no driver or
/// no devices, the registry is simply empty and matching falls back to CPU.
class CudaDeviceRegistry {
 public:
  static CudaDeviceRegistry& Instance() {
    static CudaDeviceRegistry registry;
    return registry;
  }

  int NumDevices() const { return static_cast<int>(devices_.size()); }

  int Ordinal(int idx) const { return devices_.at(idx).ordinal; }

  size_t GlobalMem(int idx) const { return devices_.at(idx).global_mem; }

 private:
  struct DeviceInfo {
    int ordinal;
    std::string name;
    size_t global_mem;
  };

  CudaDeviceRegistry() {
    int count = 0;
    // Fails harmlessly (count stays 0) when no NVIDIA driver is installed.
    if (cudaGetDeviceCount(&count) != cudaSuccess) {
      return;
    }
    for (int i = 0; i < count; ++i) {
      cudaDeviceProp prop;
      if (cudaGetDeviceProperties(&prop, i) != cudaSuccess) {
        continue;
      }
      if (prop.major * 100 + prop.minor * 10 < kMinSupportedArch) {
        std::cerr << "[CUDA] Skipping device " << prop.name
                  << " (compute capability " << prop.major << "." << prop.minor
                  << " < " << kMinSupportedArch / 100 << "."
                  << kMinSupportedArch / 10 % 10 << ")\n";
        continue;
      }
      devices_.push_back({i, prop.name, prop.totalGlobalMem});
    }
    if (!devices_.empty()) {
      std::cerr << "[CUDA] " << devices_.size()
                << " device(s) available for matching:\n";
      for (size_t i = 0; i < devices_.size(); ++i) {
        std::cerr << "  [" << i << "] " << devices_[i].name << " ("
                  << (devices_[i].global_mem >> 20) << " MB)\n";
      }
    }
  }

  std::vector<DeviceInfo> devices_;
};

// =====================================================================
// Kernels
//
// One thread per query descriptor; reference descriptors are streamed
// through shared memory in tiles that all threads of a block share —
// the same scheme as the previous OpenCL implementation.  Threads past
// n1 stay active for the cooperative tile loads and barriers but skip
// distance computation and output writes.
// =====================================================================

/// KNN2 body shared by all kernel variants.  `q` is the caller's
/// register/local copy of the query descriptor, `tile` the shared-memory
/// staging buffer for reference descriptors.
template <int NW>
__device__ void HammingKNN2Body(const uint32_t* __restrict__ q, bool valid,
                                const uint32_t* __restrict__ refs,
                                int ref_offset, int n2, uint32_t* tile,
                                int* out_idx, int* out_best, int* out_second) {
  int bd = NW * 32 + 1;
  int sd = NW * 32 + 1;
  int bi = -1;

  for (int tile_start = 0; tile_start < n2; tile_start += kHammingTileSize) {
    __syncthreads();
    const int tile_elems = kHammingTileSize * NW;
    for (int t = threadIdx.x; t < tile_elems; t += blockDim.x) {
      int row = t / NW;
      int col = t % NW;
      int j = tile_start + row;
      tile[t] =
          (j < n2) ? refs[(size_t)(ref_offset + j) * NW + col] : 0u;
    }
    __syncthreads();

    if (valid) {
      int tile_count = min(kHammingTileSize, n2 - tile_start);
      for (int t = 0; t < tile_count; ++t) {
        int dist = 0;
#pragma unroll
        for (int k = 0; k < NW; ++k) {
          dist += __popc(q[k] ^ tile[t * NW + k]);
        }
        int j = tile_start + t;
        if (dist < bd) {
          sd = bd;
          bd = dist;
          bi = j;
        } else if (dist < sd) {
          sd = dist;
        }
      }
    }
  }

  *out_idx = bi;
  *out_best = bd;
  *out_second = sd;
}

template <int NW>
__global__ void HammingKNN2Kernel(const uint32_t* __restrict__ f1,
                                  const uint32_t* __restrict__ f2,
                                  int* best_idx, int* best_dist,
                                  int* second_dist, int n1, int n2) {
  const int i = blockIdx.x * blockDim.x + threadIdx.x;
  // Do NOT return early — all threads must hit the same barriers.
  const bool valid = (i < n1);

  uint32_t q[NW];
  if (valid) {
#pragma unroll
    for (int k = 0; k < NW; ++k) {
      q[k] = f1[(size_t)i * NW + k];
    }
  }

  __shared__ uint32_t tile[kHammingTileSize * NW];

  int bi, bd, sd;
  HammingKNN2Body<NW>(q, valid, f2, 0, n2, tile, &bi, &bd, &sd);

  if (valid) {
    best_idx[i] = bi;
    best_dist[i] = bd;
    second_dist[i] = sd;
  }
}

/// Batched variant: each block processes up to kBlockSize queries of one
/// (f1, f2) pair.  block_info holds one int4 per block:
///   x = offset of this block's first query in the concatenated queries
///   y = number of valid queries in this block (≤ blockDim.x)
///   z = offset of this pair's reference block
///   w = reference count for this pair
template <int NW>
__global__ void HammingKNN2BatchedKernel(const uint32_t* __restrict__ all_f1,
                                         const uint32_t* __restrict__ all_f2,
                                         int* best_idx, int* best_dist,
                                         int* second_dist,
                                         const int4* __restrict__ block_info) {
  const int4 info = block_info[blockIdx.x];
  const bool valid = (int)threadIdx.x < info.y;
  const int out_idx = info.x + threadIdx.x;

  uint32_t q[NW];
  if (valid) {
#pragma unroll
    for (int k = 0; k < NW; ++k) {
      q[k] = all_f1[(size_t)out_idx * NW + k];
    }
  }

  __shared__ uint32_t tile[kHammingTileSize * NW];

  int bi, bd, sd;
  HammingKNN2Body<NW>(q, valid, all_f2, info.z, info.w, tile, &bi, &bd, &sd);

  if (valid) {
    best_idx[out_idx] = bi;
    best_dist[out_idx] = bd;
    second_dist[out_idx] = sd;
  }
}

/// Generic fallback for word counts without a template instantiation.
/// Uses dynamic shared memory and a fixed-capacity local query array;
/// slower than the templated kernels but correct for any n_words ≤ 32.
__global__ void HammingKNN2KernelGeneric(const uint32_t* __restrict__ f1,
                                         const uint32_t* __restrict__ f2,
                                         int* best_idx, int* best_dist,
                                         int* second_dist, int n1, int n2,
                                         int n_words) {
  extern __shared__ uint32_t tile[];

  const int i = blockIdx.x * blockDim.x + threadIdx.x;
  const bool valid = (i < n1);

  uint32_t q[kMaxWords];
  if (valid) {
    for (int k = 0; k < n_words; ++k) {
      q[k] = f1[(size_t)i * n_words + k];
    }
  }

  int bd = n_words * 32 + 1;
  int sd = n_words * 32 + 1;
  int bi = -1;

  for (int tile_start = 0; tile_start < n2; tile_start += kHammingTileSize) {
    __syncthreads();
    const int tile_elems = kHammingTileSize * n_words;
    for (int t = threadIdx.x; t < tile_elems; t += blockDim.x) {
      int row = t / n_words;
      int col = t % n_words;
      int j = tile_start + row;
      tile[t] = (j < n2) ? f2[(size_t)j * n_words + col] : 0u;
    }
    __syncthreads();

    if (valid) {
      int tile_count = min(kHammingTileSize, n2 - tile_start);
      for (int t = 0; t < tile_count; ++t) {
        int dist = 0;
        for (int k = 0; k < n_words; ++k) {
          dist += __popc(q[k] ^ tile[t * n_words + k]);
        }
        int j = tile_start + t;
        if (dist < bd) {
          sd = bd;
          bd = dist;
          bi = j;
        } else if (dist < sd) {
          sd = dist;
        }
      }
    }
  }

  if (valid) {
    best_idx[i] = bi;
    best_dist[i] = bd;
    second_dist[i] = sd;
  }
}

__global__ void HammingKNN2BatchedKernelGeneric(
    const uint32_t* __restrict__ all_f1, const uint32_t* __restrict__ all_f2,
    int* best_idx, int* best_dist, int* second_dist,
    const int4* __restrict__ block_info, int n_words) {
  extern __shared__ uint32_t tile[];

  const int4 info = block_info[blockIdx.x];
  const bool valid = (int)threadIdx.x < info.y;
  const int out_idx = info.x + threadIdx.x;

  uint32_t q[kMaxWords];
  if (valid) {
    for (int k = 0; k < n_words; ++k) {
      q[k] = all_f1[(size_t)out_idx * n_words + k];
    }
  }

  int bd = n_words * 32 + 1;
  int sd = n_words * 32 + 1;
  int bi = -1;

  for (int tile_start = 0; tile_start < info.w;
       tile_start += kHammingTileSize) {
    __syncthreads();
    const int tile_elems = kHammingTileSize * n_words;
    for (int t = threadIdx.x; t < tile_elems; t += blockDim.x) {
      int row = t / n_words;
      int col = t % n_words;
      int j = tile_start + row;
      tile[t] =
          (j < info.w) ? all_f2[(size_t)(info.z + j) * n_words + col] : 0u;
    }
    __syncthreads();

    if (valid) {
      int tile_count = min(kHammingTileSize, info.w - tile_start);
      for (int t = 0; t < tile_count; ++t) {
        int dist = 0;
        for (int k = 0; k < n_words; ++k) {
          dist += __popc(q[k] ^ tile[t * n_words + k]);
        }
        int j = tile_start + t;
        if (dist < bd) {
          sd = bd;
          bd = dist;
          bi = j;
        } else if (dist < sd) {
          sd = dist;
        }
      }
    }
  }

  if (valid) {
    best_idx[out_idx] = bi;
    best_dist[out_idx] = bd;
    second_dist[out_idx] = sd;
  }
}

// =====================================================================
// Host-side launch helpers
// =====================================================================

void ValidateWordCount(int n_words) {
  if (n_words <= 0 || n_words > kMaxWords) {
    throw std::runtime_error("Unsupported binary descriptor word count: " +
                             std::to_string(n_words) + " (max " +
                             std::to_string(kMaxWords) + ")");
  }
}

/// Select the device backing registry index `device_idx` for this thread.
void ActivateDevice(int device_idx) {
  auto& registry = CudaDeviceRegistry::Instance();
  if (device_idx < 0 || device_idx >= registry.NumDevices()) {
    throw std::runtime_error("Invalid CUDA device index");
  }
  CUDA_CHECK(cudaSetDevice(registry.Ordinal(device_idx)));
}

void LaunchKNN2(const uint32_t* d_f1, const uint32_t* d_f2, int* d_best_idx,
                int* d_best_dist, int* d_second_dist, int n1, int n2,
                int n_words) {
  const int blocks = (n1 + kBlockSize - 1) / kBlockSize;
  switch (n_words) {
    case 2:
      HammingKNN2Kernel<2><<<blocks, kBlockSize>>>(
          d_f1, d_f2, d_best_idx, d_best_dist, d_second_dist, n1, n2);
      break;
    case 4:
      HammingKNN2Kernel<4><<<blocks, kBlockSize>>>(
          d_f1, d_f2, d_best_idx, d_best_dist, d_second_dist, n1, n2);
      break;
    case 8:
      HammingKNN2Kernel<8><<<blocks, kBlockSize>>>(
          d_f1, d_f2, d_best_idx, d_best_dist, d_second_dist, n1, n2);
      break;
    default: {
      const size_t shared_bytes =
          (size_t)kHammingTileSize * n_words * sizeof(uint32_t);
      HammingKNN2KernelGeneric<<<blocks, kBlockSize, shared_bytes>>>(
          d_f1, d_f2, d_best_idx, d_best_dist, d_second_dist, n1, n2, n_words);
      break;
    }
  }
  CUDA_CHECK(cudaGetLastError());
  CUDA_CHECK(cudaDeviceSynchronize());
}

void LaunchKNN2Batched(const uint32_t* d_f1, const uint32_t* d_f2,
                       int* d_best_idx, int* d_best_dist, int* d_second_dist,
                       const int4* d_block_info, int total_blocks,
                       int n_words) {
  switch (n_words) {
    case 2:
      HammingKNN2BatchedKernel<2><<<total_blocks, kBlockSize>>>(
          d_f1, d_f2, d_best_idx, d_best_dist, d_second_dist, d_block_info);
      break;
    case 4:
      HammingKNN2BatchedKernel<4><<<total_blocks, kBlockSize>>>(
          d_f1, d_f2, d_best_idx, d_best_dist, d_second_dist, d_block_info);
      break;
    case 8:
      HammingKNN2BatchedKernel<8><<<total_blocks, kBlockSize>>>(
          d_f1, d_f2, d_best_idx, d_best_dist, d_second_dist, d_block_info);
      break;
    default: {
      const size_t shared_bytes =
          (size_t)kHammingTileSize * n_words * sizeof(uint32_t);
      HammingKNN2BatchedKernelGeneric<<<total_blocks, kBlockSize,
                                        shared_bytes>>>(
          d_f1, d_f2, d_best_idx, d_best_dist, d_second_dist, d_block_info,
          n_words);
      break;
    }
  }
  CUDA_CHECK(cudaGetLastError());
  CUDA_CHECK(cudaDeviceSynchronize());
}

}  // namespace

// =====================================================================
// Public interface (matching_gpu_cuda.h)
// =====================================================================

bool CudaAvailable() {
  return CudaDeviceRegistry::Instance().NumDevices() > 0;
}

int CudaNumDevices() { return CudaDeviceRegistry::Instance().NumDevices(); }

size_t CudaDeviceGlobalMem(int device_idx) {
  return CudaDeviceRegistry::Instance().GlobalMem(device_idx);
}

KNN2ResultInt HammingKNN2(int device_idx, const uint32_t* f1, int n1,
                          const uint32_t* f2, int n2, int n_words) {
  ValidateWordCount(n_words);
  ActivateDevice(device_idx);

  const size_t f1_bytes = (size_t)n1 * n_words * sizeof(uint32_t);
  const size_t f2_bytes = (size_t)n2 * n_words * sizeof(uint32_t);

  DeviceBuffer d_f1(f1_bytes);
  DeviceBuffer d_f2(f2_bytes);
  DeviceBuffer d_best_idx((size_t)n1 * sizeof(int));
  DeviceBuffer d_best_dist((size_t)n1 * sizeof(int));
  DeviceBuffer d_second_dist((size_t)n1 * sizeof(int));

  CUDA_CHECK(cudaMemcpy(d_f1.as<uint32_t>(), f1, f1_bytes,
                        cudaMemcpyHostToDevice));
  CUDA_CHECK(cudaMemcpy(d_f2.as<uint32_t>(), f2, f2_bytes,
                        cudaMemcpyHostToDevice));

  LaunchKNN2(d_f1.as<uint32_t>(), d_f2.as<uint32_t>(), d_best_idx.as<int>(),
             d_best_dist.as<int>(), d_second_dist.as<int>(), n1, n2, n_words);

  KNN2ResultInt result;
  result.best_idx.resize(n1);
  result.best_dist.resize(n1);
  result.second_dist.resize(n1);
  CUDA_CHECK(cudaMemcpy(result.best_idx.data(), d_best_idx.as<int>(),
                        (size_t)n1 * sizeof(int), cudaMemcpyDeviceToHost));
  CUDA_CHECK(cudaMemcpy(result.best_dist.data(), d_best_dist.as<int>(),
                        (size_t)n1 * sizeof(int), cudaMemcpyDeviceToHost));
  CUDA_CHECK(cudaMemcpy(result.second_dist.data(), d_second_dist.as<int>(),
                        (size_t)n1 * sizeof(int), cudaMemcpyDeviceToHost));
  return result;
}

KNN2ResultInt HammingKNN2Batched(int device_idx, const uint32_t* queries,
                                 int total_queries, const uint32_t* refs,
                                 int total_refs,
                                 const std::vector<PairOffsets>& pairs,
                                 int n_words) {
  ValidateWordCount(n_words);
  ActivateDevice(device_idx);

  // One int4 per block: {first query offset, valid queries, ref offset,
  // ref count}.  Same partitioning as the previous OpenCL wg_info scheme,
  // with a fixed block size.
  std::vector<int4> block_info;
  for (const auto& p : pairs) {
    for (int off = 0; off < p.n1; off += kBlockSize) {
      int4 info;
      info.x = p.f1_offset + off;
      info.y = std::min(kBlockSize, p.n1 - off);
      info.z = p.f2_offset;
      info.w = p.n2;
      block_info.push_back(info);
    }
  }

  KNN2ResultInt result;
  result.best_idx.assign(total_queries, -1);
  result.best_dist.assign(total_queries, n_words * 32 + 1);
  result.second_dist.assign(total_queries, n_words * 32 + 1);
  if (block_info.empty() || total_queries == 0) {
    return result;
  }

  const size_t q_bytes = (size_t)total_queries * n_words * sizeof(uint32_t);
  const size_t r_bytes = (size_t)total_refs * n_words * sizeof(uint32_t);
  const size_t info_bytes = block_info.size() * sizeof(int4);

  DeviceBuffer d_q(q_bytes);
  DeviceBuffer d_r(std::max(r_bytes, (size_t)1));
  DeviceBuffer d_best_idx((size_t)total_queries * sizeof(int));
  DeviceBuffer d_best_dist((size_t)total_queries * sizeof(int));
  DeviceBuffer d_second_dist((size_t)total_queries * sizeof(int));
  DeviceBuffer d_block_info(info_bytes);

  CUDA_CHECK(cudaMemcpy(d_q.as<uint32_t>(), queries, q_bytes,
                        cudaMemcpyHostToDevice));
  if (r_bytes > 0) {
    CUDA_CHECK(cudaMemcpy(d_r.as<uint32_t>(), refs, r_bytes,
                          cudaMemcpyHostToDevice));
  }
  CUDA_CHECK(cudaMemcpy(d_block_info.as<int4>(), block_info.data(), info_bytes,
                        cudaMemcpyHostToDevice));

  // Queries not covered by any block (pairs with n1 == 0 contribute no
  // blocks) keep their host-side initial values; blocks with n2 == 0
  // write bi = -1, so every query row ends up well-defined.
  CUDA_CHECK(cudaMemcpy(d_best_idx.as<int>(), result.best_idx.data(),
                        (size_t)total_queries * sizeof(int),
                        cudaMemcpyHostToDevice));

  LaunchKNN2Batched(d_q.as<uint32_t>(), d_r.as<uint32_t>(),
                    d_best_idx.as<int>(), d_best_dist.as<int>(),
                    d_second_dist.as<int>(), d_block_info.as<int4>(),
                    (int)block_info.size(), n_words);

  CUDA_CHECK(cudaMemcpy(result.best_idx.data(), d_best_idx.as<int>(),
                        (size_t)total_queries * sizeof(int),
                        cudaMemcpyDeviceToHost));
  CUDA_CHECK(cudaMemcpy(result.best_dist.data(), d_best_dist.as<int>(),
                        (size_t)total_queries * sizeof(int),
                        cudaMemcpyDeviceToHost));
  CUDA_CHECK(cudaMemcpy(result.second_dist.data(), d_second_dist.as<int>(),
                        (size_t)total_queries * sizeof(int),
                        cudaMemcpyDeviceToHost));
  return result;
}

}  // namespace gpu
}  // namespace features
