#pragma once
#include <cstdint>
#include <cstdlib>
#include <malloc.h>

#include "bridge.h"

namespace gbbs {
class BitArray {
private:
  uint64_t len;
  uint32_t *array;

  static inline uint64_t bit_array_size(uint64_t size) {
    if (size == 0) {
      return 0;
    }
    if (size < 32) {
      size = 32;
    }
    uint64_t n = size / 32;
    if (n * 32 < size) {
      n += 1;
    }
    return n * 4;
  }

public:
  //   BitArray(uint32_t *arr, uint64_t size)
  //       : array(arr), len(size), to_free(false) {}
  BitArray() : len(0), array(nullptr) {}
  explicit BitArray(uint64_t size) {
    uint64_t n = bit_array_size(size);
    array = (uint32_t *)memalign(32, n);
    len = n * 8;
    parallel_for(0, len / 32, [&](size_t i) { array[i] = 0; });
  }
  explicit BitArray(uint64_t size, bool set_all) {
    uint64_t n = bit_array_size(size);
    array = (uint32_t *)memalign(32, n);
    len = n * 8;
    if (set_all) {
      parallel_for(0, len / 32, [&](size_t i) { array[i] = 0xFFFFFFFFU; });
    } else {
      parallel_for(0, len / 32, [&](size_t i) { array[i] = 0; });
    }
  }
  static BitArray uninitialized(uint64_t size) {
    BitArray ba;
    uint64_t n = bit_array_size(size);
    ba.array = (uint32_t *)memalign(32, n);
    ba.len = n * 8;
    return ba;
  }

  BitArray(const BitArray &other) : len(other.len) {
    array = (uint32_t *)memalign(32, len / 8);
    parallel_for(0, len / 32, [&](size_t i) { array[i] = other.array[i]; });
  }

  BitArray &operator=(BitArray &&other) {
    if (this != &other) {
      free(array);
      len = other.len;
      array = other.array;
      other.array = nullptr;
      other.len = 0;
    }
    return *this;
  }

  BitArray(BitArray &&other) : len(other.len), array(other.array) {
    other.array = nullptr;
    other.len = 0;
  }
  void clear() const {
    parallel_for(0, len / 32, [&](size_t i) { array[i] = 0; });
  }

  ~BitArray() { free(array); }
  [[nodiscard]] bool get(uint64_t i) const {
    return (array[i / 32] >> i % 32) & 1U;
  }
  [[nodiscard]] bool operator[](size_t i) const { return get(i); }
  void set(uint64_t i) const { array[i / 32] |= (1U << i % 32); }
  void set_atomic(uint64_t i) const {
    __atomic_fetch_or(&array[i / 32], 1U << (i % 32), __ATOMIC_RELAXED);
  }
  void flip(uint64_t i) const { array[i / 32] ^= (1U << i % 32); }

  [[nodiscard]] uint64_t count() const {
    uint64_t count = 0;
    for (uint64_t i = 0; i < len / 32; i++) {
      count += __builtin_popcount(array[i]);
    }
    return count;
  }
  [[nodiscard]] bool non_empty() const {
    for (uint64_t i = 0; i < len / 32; i++) {
      if (array[i] > 0) {
        return true;
      }
    }
    return false;
  }
  template <bool parallel, class F> void map(F &f) const {
    if constexpr (parallel) {
      aligned_parallel_for(0, len, 256, [&](size_t i) {
        if (get(i)) {
          f(i);
        }
      });
    } else {
      for (uint64_t i = 0; i < len; i++) {
        if (get(i)) {
          f(i);
        }
      }
    }
  }
  size_t size() const { return len; }
};
/*
// used the check the impact of the bit array
class BitArray {
private:
  uint64_t len;
  bool *array;

  static inline uint64_t bit_array_size(uint64_t size) {
    if (size == 0) {
      return 0;
    }
    if (size < 32) {
      size = 32;
    }
    uint64_t n = size / 32;
    if (n * 32 < size) {
      n += 1;
    }
    return n * 32;
  }

public:
  //   BitArray(uint32_t *arr, uint64_t size)
  //       : array(arr), len(size), to_free(false) {}
  BitArray() : len(0), array(nullptr) {}
  explicit BitArray(uint64_t size) {
    uint64_t n = bit_array_size(size);
    array = (bool *)memalign(32, n);
    len = size;
    parallel_for(0, len, [&](size_t i) { array[i] = 0; });
  }
  explicit BitArray(uint64_t size, bool set_all) {
    uint64_t n = bit_array_size(size);
    array = (bool *)memalign(32, n);
    len = size;
    if (set_all) {
      parallel_for(0, len, [&](size_t i) { array[i] = true; });
    } else {
      parallel_for(0, len, [&](size_t i) { array[i] = false; });
    }
  }
  static BitArray uninitialized(uint64_t size) {
    BitArray ba;
    uint64_t n = bit_array_size(size);
    ba.array = (bool *)memalign(32, n);
    ba.len = size;
    return ba;
  }

  BitArray(const BitArray &other) : len(other.len) {
    array = (bool *)memalign(32, bit_array_size(len));
    parallel_for(0, len, [&](size_t i) { array[i] = other.array[i]; });
  }

  BitArray &operator=(BitArray &&other) {
    if (this != &other) {
      free(array);
      len = other.len;
      array = other.array;
      other.array = nullptr;
      other.len = 0;
    }
    return *this;
  }

  BitArray(BitArray &&other) : len(other.len), array(other.array) {
    other.array = nullptr;
    other.len = 0;
  }
  void clear() const {
    parallel_for(0, len, [&](size_t i) { array[i] = 0; });
  }

  ~BitArray() { free(array); }
  [[nodiscard]] bool get(uint64_t i) const {
    return array[i];
  }
  [[nodiscard]] bool operator[](size_t i) const { return get(i); }
  void set(uint64_t i) const { array[i] = true; }
  void set_atomic(uint64_t i) const {
    array[i] = true;
  }
  void flip(uint64_t i) const { array[i] = !array[i]; }

  [[nodiscard]] uint64_t count() const {
    uint64_t count = 0;
    for (uint64_t i = 0; i < len; i++) {
      count += array[i];
    }
    return count;
  }
  [[nodiscard]] bool non_empty() const {
    for (uint64_t i = 0; i < len; i++) {
      if (array[i] > 0) {
        return true;
      }
    }
    return false;
  }
  template <bool parallel, class F> void map(F &f) const {
    if constexpr (parallel) {
      aligned_parallel_for(0, len, 256, [&](size_t i) {
        if (get(i)) {
          f(i);
        }
      });
    } else {
      for (uint64_t i = 0; i < len; i++) {
        if (get(i)) {
          f(i);
        }
      }
    }
  }
  size_t size() const { return len; }
};
*/
} // namespace gbbs