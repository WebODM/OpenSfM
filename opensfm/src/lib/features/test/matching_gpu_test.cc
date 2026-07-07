#include <features/matching_gpu.h>
#include <foundation/python_types.h>
#include <gtest/gtest.h>
#include <pybind11/embed.h>
#include <pybind11/numpy.h>

#include <cstring>
#include <random>
#include <set>
#include <utility>
#include <vector>

namespace py = pybind11;

namespace {

// The GPU matcher API takes/returns numpy arrays, so the embedded Python
// interpreter must be initialized before any py::array_t is constructed.
py::scoped_interpreter kPythonInterpreter{};

constexpr int kWords = 4;  // 128-bit binary descriptors

// Generate N random packed binary descriptors of kWords uint32 words.
std::vector<uint32_t> RandomBinaryDescriptors(int n, unsigned seed) {
  std::mt19937 rng(seed);
  std::vector<uint32_t> data(n * kWords);
  for (auto& v : data) {
    v = rng();
  }
  return data;
}

// Copy descriptors f1[from..from+n) into f2[to..to+n) with a single bit
// flipped, planting near-duplicate matches.
void PlantMatches(const std::vector<uint32_t>& f1, std::vector<uint32_t>& f2,
                  int n, int from, int to) {
  for (int i = 0; i < n; ++i) {
    std::memcpy(&f2[(to + i) * kWords], &f1[(from + i) * kWords],
                kWords * sizeof(uint32_t));
    f2[(to + i) * kWords] ^= 1u;  // Hamming distance 1
  }
}

int PopCount(uint32_t v) {
  int c = 0;
  for (; v != 0; v &= v - 1) {
    ++c;
  }
  return c;
}

struct Match {
  int query;
  int ref;
};

// CPU brute-force Hamming 2-NN with Lowe's ratio test, for reference.
std::vector<Match> CpuHamming(const uint32_t* f1, int n1, const uint32_t* f2,
                              int n2, float lowes_ratio) {
  std::vector<Match> matches;
  for (int i = 0; i < n1; ++i) {
    int bd = kWords * 32 + 1, sd = kWords * 32 + 1, bi = -1;
    for (int j = 0; j < n2; ++j) {
      int dist = 0;
      for (int k = 0; k < kWords; ++k) {
        dist += PopCount(f1[i * kWords + k] ^ f2[j * kWords + k]);
      }
      if (dist < bd) {
        sd = bd;
        bd = dist;
        bi = j;
      } else if (dist < sd) {
        sd = dist;
      }
    }
    if (bi >= 0 && (float)bd < lowes_ratio * (float)sd) {
      matches.push_back({i, bi});
    }
  }
  return matches;
}

py::array_t<uint32_t> ToArray(const std::vector<uint32_t>& data, int n) {
  py::array_t<uint32_t> arr({n, kWords});
  if (!data.empty()) {
    std::memcpy(arr.mutable_data(), data.data(),
                data.size() * sizeof(uint32_t));
  }
  return arr;
}

std::set<std::pair<int, int>> ToSet(const py::array_t<int>& matches) {
  std::set<std::pair<int, int>> result;
  auto r = matches.unchecked<2>();
  for (py::ssize_t k = 0; k < matches.shape(0); ++k) {
    result.emplace(r(k, 0), r(k, 1));
  }
  return result;
}

}  // namespace

TEST(MatchingGPU, Availability) {
  // Just check that the function doesn't crash.
  bool avail = features::gpu_matching_available();
  if (!avail) {
    GTEST_SUCCEED() << "CUDA not available, skipping remaining tests";
    return;
  }
  EXPECT_GE(features::gpu_num_devices(), 1);
}

TEST(MatchingGPU, BasicHammingMatching) {
  if (!features::gpu_matching_available()) {
    GTEST_SUCCEED() << "CUDA not available";
    return;
  }

  constexpr int n1 = 200;
  constexpr int n2 = 300;
  constexpr float ratio = 0.8f;

  auto data1 = RandomBinaryDescriptors(n1, 42);
  auto data2 = RandomBinaryDescriptors(n2, 123);
  PlantMatches(data1, data2, 10, 0, 50);

  auto result =
      features::match_hamming_gpu(ToArray(data1, n1), ToArray(data2, n2),
                                  ratio);
  ASSERT_GE(result.shape(0), 10);
  ASSERT_EQ(result.shape(1), 2);

  auto matches = ToSet(result);
  for (int i = 0; i < 10; ++i) {
    EXPECT_TRUE(matches.count({i, 50 + i}));
  }
}

TEST(MatchingGPU, SymmetricHammingMatching) {
  if (!features::gpu_matching_available()) {
    GTEST_SUCCEED() << "CUDA not available";
    return;
  }

  constexpr int n1 = 150;
  constexpr int n2 = 200;
  constexpr float ratio = 0.8f;

  auto data1 = RandomBinaryDescriptors(n1, 77);
  auto data2 = RandomBinaryDescriptors(n2, 88);
  PlantMatches(data1, data2, 5, 10, 30);

  auto f1 = ToArray(data1, n1);
  auto f2 = ToArray(data2, n2);

  auto sym = features::match_hamming_gpu_symmetric(f1, f2, ratio);
  ASSERT_GE(sym.shape(0), 5);
  ASSERT_EQ(sym.shape(1), 2);

  auto sym_set = ToSet(sym);
  for (int i = 0; i < 5; ++i) {
    EXPECT_TRUE(sym_set.count({10 + i, 30 + i}));
  }

  // Symmetric should return fewer or equal matches than asymmetric.
  auto asym = features::match_hamming_gpu(f1, f2, ratio);
  EXPECT_LE(sym.shape(0), asym.shape(0));
}

TEST(MatchingGPU, EmptyInput) {
  if (!features::gpu_matching_available()) {
    GTEST_SUCCEED() << "CUDA not available";
    return;
  }

  py::array_t<uint32_t> f1(std::vector<int>{0, kWords});
  auto data2 = RandomBinaryDescriptors(100, 1);
  auto f2 = ToArray(data2, 100);

  auto result = features::match_hamming_gpu(f1, f2, 0.8f);
  EXPECT_EQ(result.shape(0), 0);
  EXPECT_EQ(result.shape(1), 2);

  auto result2 = features::match_hamming_gpu_symmetric(f1, f2, 0.8f);
  EXPECT_EQ(result2.shape(0), 0);
}

TEST(MatchingGPU, ConsistencyWithCPU) {
  if (!features::gpu_matching_available()) {
    GTEST_SUCCEED() << "CUDA not available";
    return;
  }

  constexpr int n1 = 100;
  constexpr int n2 = 120;
  constexpr float ratio = 0.75f;

  auto data1 = RandomBinaryDescriptors(n1, 1234);
  auto data2 = RandomBinaryDescriptors(n2, 5678);

  // Plant exact copies to ensure matches.
  for (int i = 0; i < 8; ++i) {
    std::memcpy(&data2[(i + 10) * kWords], &data1[i * kWords],
                kWords * sizeof(uint32_t));
  }

  auto gpu_result = features::match_hamming_gpu(ToArray(data1, n1),
                                                ToArray(data2, n2), ratio);
  auto cpu_matches = CpuHamming(data1.data(), n1, data2.data(), n2, ratio);

  // GPU and CPU brute-force should produce identical results.
  ASSERT_EQ(gpu_result.shape(0), (py::ssize_t)cpu_matches.size());

  auto r = gpu_result.unchecked<2>();
  for (size_t k = 0; k < cpu_matches.size(); ++k) {
    EXPECT_EQ(r(k, 0), cpu_matches[k].query);
    EXPECT_EQ(r(k, 1), cpu_matches[k].ref);
  }
}

TEST(MatchingGPU, BatchedMatchesSinglePair) {
  if (!features::gpu_matching_available()) {
    GTEST_SUCCEED() << "CUDA not available";
    return;
  }

  constexpr float ratio = 0.8f;
  const int sizes1[] = {120, 0, 350, 64};
  const int sizes2[] = {200, 80, 0, 512};

  py::list f1_list, f2_list;
  std::vector<py::array_t<uint32_t>> arrays1, arrays2;
  for (int p = 0; p < 4; ++p) {
    auto d1 = RandomBinaryDescriptors(sizes1[p], 100 + p);
    auto d2 = RandomBinaryDescriptors(sizes2[p], 200 + p);
    if (sizes1[p] >= 20 && sizes2[p] >= 40) {
      PlantMatches(d1, d2, 10, 0, 20);
    }
    arrays1.push_back(ToArray(d1, sizes1[p]));
    arrays2.push_back(ToArray(d2, sizes2[p]));
    f1_list.append(arrays1.back());
    f2_list.append(arrays2.back());
  }

  auto batch_results =
      features::match_hamming_gpu_batch_symmetric(f1_list, f2_list, ratio);
  ASSERT_EQ(py::len(batch_results), 4);

  // Each pair's batched result must equal the single-pair symmetric result.
  for (int p = 0; p < 4; ++p) {
    auto batched = batch_results[p].cast<py::array_t<int>>();
    if (sizes1[p] == 0 || sizes2[p] == 0) {
      EXPECT_EQ(batched.shape(0), 0) << "pair " << p;
      continue;
    }
    auto single =
        features::match_hamming_gpu_symmetric(arrays1[p], arrays2[p], ratio);
    EXPECT_EQ(ToSet(batched), ToSet(single)) << "pair " << p;
  }
}
