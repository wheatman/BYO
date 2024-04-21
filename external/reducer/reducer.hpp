#pragma once

#include "parlay/parallel.h"
#include "parlay/primitives.h"
#include "parlay/sequence.h"

#include <algorithm>
#include <cstring>
#include <functional>
#include <limits>
#include <type_traits>
#include <vector>

#if CILK == 1
#include <cilk/cilksan.h>
#endif

template <class F> class Reducer {

#ifdef __cpp_lib_hardware_interference_size
  using std::hardware_constructive_interference_size;
  using std::hardware_destructive_interference_size;
#else
  // 64 bytes on x86-64 │ L1_CACHE_BYTES │ L1_CACHE_SHIFT │ __cacheline_aligned
  // │
  // ...
  static constexpr std::size_t hardware_constructive_interference_size = 64;
  static constexpr std::size_t hardware_destructive_interference_size = 64;
#endif

  struct aligned_f {
    alignas(hardware_destructive_interference_size) F f;
    aligned_f(F f_) : f(f_) {}
  };
  std::vector<aligned_f> data;

#if CILK == 1
  // so cilksan doesn't report races on accesses to the vector which I make sure
  // are fine by using getWorkerNum()
  Cilksan_fake_mutex fake_lock;
#endif

public:
  Reducer(F identity = {}) { data.resize(parlay::num_workers(), identity); }
  void update(F new_values) {
    int worker_num = parlay::worker_id();
#if CILK == 1
    Cilksan_fake_lock_guard guad(&fake_lock);
#endif
    data[worker_num].f.update(new_values);
  }
  F get() const {
    F output;
    for (const auto &d : data) {
      output.update(d.f);
    }
    return output;
  }
};

template <class Monoid> class Reducer_with_object {

#ifdef __cpp_lib_hardware_interference_size
  using std::hardware_constructive_interference_size;
  using std::hardware_destructive_interference_size;
#else
  // 64 bytes on x86-64 │ L1_CACHE_BYTES │ L1_CACHE_SHIFT │ __cacheline_aligned
  // │
  // ...
  static constexpr std::size_t hardware_constructive_interference_size = 64;
  static constexpr std::size_t hardware_destructive_interference_size = 64;
#endif

  using F = typename Monoid::T;

  struct aligned_f {
    alignas(hardware_destructive_interference_size) F f;
    aligned_f(F f_) : f(f_) {}
  };
  std::vector<aligned_f> data;
  Monoid monoid;

#if CILK == 1
  // so cilksan doesn't report races on accesses to the vector which I make sure
  // are fine by using getWorkerNum()
  Cilksan_fake_mutex fake_lock;
#endif

public:
  Reducer_with_object(Monoid m) : monoid(m) {
    data.resize(parlay::num_workers(), m.identity);
  }
  void update(F new_value) {
    int worker_num = parlay::worker_id();
#if CILK == 1
    Cilksan_fake_lock_guard guad(&fake_lock);
#endif
    data[worker_num].f = monoid(data[worker_num].f, new_value);
  }
  F get() const {
    F output = monoid.identity;
    for (const auto &d : data) {
      output = monoid(output, d.f);
    }
    return output;
  }
};

template <class T> class Reducer_sum {
  static_assert(std::is_integral<T>::value, "Integral required.");
  struct F {
    T value = 0;
    void update(const F new_value) { value += new_value.value; }
    F(T t) : value(t) {}
    F() {}
  };
  Reducer<F> reducer;

public:
  Reducer_sum(T initial_value = {}) { add(initial_value); }
  void add(T new_value) { reducer.update(new_value); }
  void inc() { reducer.update(1); }
  T get() const { return reducer.get().value; }
  Reducer_sum &operator++() {
    inc();
    return *this;
  }
  Reducer_sum &operator--() {
    add(-1);
    return *this;
  }

  Reducer_sum &operator-=(T new_value) {
    add(-new_value);
    return *this;
  }
  Reducer_sum &operator+=(T new_value) {
    add(+new_value);
    return *this;
  }

  friend bool operator==(const Reducer_sum &lhs, const Reducer_sum &rhs) {
    return lhs.get() == rhs.get();
  }
  operator T() const { return get(); }
};

template <class T> class Reducer_max {
  struct F {
    T value = std::numeric_limits<T>::min();
    void update(const F new_value) { value = std::max(value, new_value.value); }
    F(T t) : value(t) {}
    F() {}
  };
  Reducer<F> reducer;

public:
  Reducer_max() {}
  void update(T new_value) { reducer.update(new_value); }
  T get() const { return reducer.get().value; }
};

template <class T> class Reducer_Vector {

#ifdef __cpp_lib_hardware_interference_size
  using std::hardware_constructive_interference_size;
  using std::hardware_destructive_interference_size;
#else
  // 64 bytes on x86-64 │ L1_CACHE_BYTES │ L1_CACHE_SHIFT │ __cacheline_aligned
  // │
  // ...
  static constexpr std::size_t hardware_constructive_interference_size = 128;
  static constexpr std::size_t hardware_destructive_interference_size = 128;
#endif

  struct aligned_f {
    alignas(hardware_destructive_interference_size) std::vector<T> f;
  };
  std::vector<aligned_f> data;

public:
  using elements_type = T;
  Reducer_Vector() { data.resize(parlay::num_workers()); }

  Reducer_Vector(std::vector<T> &start) {
    data.resize(parlay::worker_id());
    data[0].f = std::move(start);
  }

  template <typename F> void push_back(F arg) {
    static_assert(std::is_constructible_v<T, F>);
    int worker_num = parlay::worker_id();
    data[worker_num].f.emplace_back(arg);
  }
  void push_back(T arg) {
    int worker_num = parlay::worker_id();
    data[worker_num].f.push_back(arg);
  }
  std::vector<T> get_sorted() const {
    std::vector<size_t> lengths(data.size() + 1);
    for (size_t i = 1; i <= data.size(); i++) {
      lengths[i] += lengths[i - 1] + data[i - 1].f.size();
    }
    std::vector<T> output(lengths[data.size()]);
    parallel_for(0, data.size(), [&](size_t i) {
      std::memcpy(output.data() + lengths[i], data[i].f.data(),
                  data[i].f.size() * sizeof(T));
    }, 1);
    parlay::sort(output.begin());
    return output;
  }

  parlay::sequence<T> get_sorted_sequence() const {
    std::vector<size_t> lengths(data.size() + 1);
    for (size_t i = 1; i <= data.size(); i++) {
      lengths[i] += lengths[i - 1] + data[i - 1].f.size();
    }
    parlay::sequence<T> output =
        parlay::sequence<T>::uninitialized(lengths[data.size()]);
    parlay::parallel_for(0, data.size(), [&](size_t i) {
      std::uninitialized_move(data[i].f.begin(), data[i].f.end(),
                              output.begin() + lengths[i]);
    }, 1);
    parlay::sort(output);
    return output;
  }

  parlay::sequence<T> get_sequence() const {
    std::vector<size_t> lengths(data.size() + 1);
    for (size_t i = 1; i <= data.size(); i++) {
      lengths[i] += lengths[i - 1] + data[i - 1].f.size();
    }
    parlay::sequence<T> output =
        parlay::sequence<T>::uninitialized(lengths[data.size()]);
    parlay::parallel_for(0, data.size(), [&](size_t i) {
      std::uninitialized_move(data[i].f.begin(), data[i].f.end(),
                              output.begin() + lengths[i]);
    }, 1);
    return output;
  }

  std::vector<T> get() const {
    if (data.size() == 0) {
      return {};
    }
    std::vector<size_t> lengths(data.size() + 1);
    for (size_t i = 1; i <= data.size(); i++) {
      lengths[i] += lengths[i - 1] + data[i - 1].f.size();
    }
    if (lengths[data.size()] == 0) {
      return {};
    }
    std::vector<T> output(lengths[data.size()]);
    parlay::parallel_for(0, data.size(), [&](size_t i) {
      std::memcpy(output.data() + lengths[i], data[i].f.data(),
                  data[i].f.size() * sizeof(T));
    }, 1);
    return output;
  }

  template <typename F> void for_each(F f) const {
    parlay::parallel_for(0, data.size(), [&](size_t i) {
      parlay::parallel_for(0, data[i].f.size(),
                           [&](size_t j) { f(data[i].f[j]); });
    });
  }
  template <typename F> void serial_for_each(F f) const {
    for (auto &d : data) {
      for (auto &e : d.f) {
        f(e);
      }
    }
  }
  template <typename C, typename R, typename K>
  K find_first_match(C c, R r, K default_return) {
    for (auto &d : data) {
      for (auto &e : d.f) {
        if (c(e)) {
          return r(e);
        }
      }
    }
    return default_return;
  }

  size_t size() const {
    size_t n = 0;
    for (auto &d : data) {
      n += d.f.size();
    }
    return n;
  }
  bool empty() const {
    for (auto &d : data) {
      if (!d.f.empty()) {
        return false;
      }
    }
    return true;
  }

  // the following parameters can be tuned
  static constexpr const size_t _cs_seq_threshold = 2048;
  static constexpr const size_t _cs_max_blocks = 512;

  // Sequential base case for the no-transpose count sort.
  template <typename b_size_t, typename s_size_t, typename E, typename I,
            typename F>
  inline void _seq_count_sort(I *In, E *Out, F &get_key, s_size_t region_length,
                              s_size_t *counts, s_size_t num_buckets) {
    auto offsets = parlay::sequence<s_size_t>::uninitialized(num_buckets);
    auto tmp = parlay::sequence<b_size_t>::uninitialized(region_length);

    for (s_size_t i = 0; i < num_buckets; i++) {
      offsets[i] = 0;
    }
    for (s_size_t j = 0; j < region_length; j++) {
      s_size_t k = tmp[j] = get_key(In[j]);
      offsets[k]++;
    }
    size_t s = 0;
    for (s_size_t i = 0; i < num_buckets; i++) {
      s += offsets[i];
      offsets[i] = s;
      counts[i] = s;
    }
    for (long j = ((long)region_length) - 1; j >= 0; j--) {
      s_size_t k = --offsets[tmp[j]];
      // needed for types with self defined assignment or initialization
      // otherwise equivalent to: Out[k+start] = In[j+start];
      parlay::assign_uninitialized(Out[k], In[j]);
    }
  }

  template <typename b_size_t, typename s_size_t, typename E, typename F>
  inline std::tuple<parlay::sequence<E>, parlay::sequence<s_size_t>, s_size_t,
                    s_size_t>
  _count_sort(F &get_key, s_size_t num_buckets) {
    size_t num_blocks = data.size();
    std::vector<size_t> lengths(data.size() + 1);
    for (size_t i = 1; i <= data.size(); i++) {
      lengths[i] += lengths[i - 1] + data[i - 1].f.size();
    }
    size_t n = lengths[data.size()];

    if (n < _cs_seq_threshold) {
      auto counts = parlay::sequence<s_size_t>::uninitialized(num_buckets + 1);
      auto seq = get_sequence();
      auto B = parlay::sequence<E>::uninitialized(n);
      _seq_count_sort<b_size_t>(seq.begin(), B.begin(), get_key, n,
                                counts.begin(), num_buckets);
      return std::make_tuple(B, counts, (s_size_t)1, num_buckets + 1);
    }

    s_size_t m = num_blocks * num_buckets;

    auto B = parlay::sequence<E>::uninitialized(n);
    auto counts = parlay::sequence<s_size_t>::uninitialized(m);

    parlay::parallel_for(0, data.size(), [&](size_t i) {
      _seq_count_sort<b_size_t>(data[i].f.data(), B.begin() + lengths[i],
                                get_key, data[i].f.size(),
                                counts.begin() + i * num_buckets, num_buckets);
    }, 1);
    return std::make_tuple(std::move(B), std::move(counts), num_blocks, m);
  }

  inline size_t block_size(size_t block_number) const {
    return data[block_number].f.size();
  }

  std::vector<size_t> prefix_sum_array() const {
    std::vector<size_t> lengths(data.size() + 1);
    for (size_t i = 1; i <= data.size(); i++) {
      lengths[i] += lengths[i - 1] + data[i - 1].f.size();
    }
    return lengths;
  }

  T get_by_index(size_t index) const {
    assert(index < size());
    for (size_t i = 0; i < data.size(); i++) {
      if (index < data[i].f.size()) {
        return data[i].f[index];
      } else {
        index -= data[i].f.size();
      }
    }
    // should never happen
    return 0;
  }
};
