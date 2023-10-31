#pragma once

#include <vector>

#include "gbbs/gbbs.h"
#include "reducer.hpp"

namespace gbbs {

template <typename T, size_t alignment> class wide_int {
  alignas(alignment) T value;

public:
  wide_int(T v) : value(v) {}
  wide_int() : value(0) {}
  operator T() const { return value; }

  auto operator<=>(const wide_int &b) const { return value <=> b.value; }
};

template <class Graph, class T>
size_t Test_Map_Remote(const Graph &G, const std::vector<T> &remote_data) {
  Reducer_sum<size_t> count;

  parallel_for(0, G.N(), [&count, &G, &remote_data](size_t i) {
    size_t local_count = 0;
    G.map_neighbors(i, [&](auto src, auto dest, auto &val) {
      local_count += remote_data[dest];
    });
    count += local_count;
  });
  return count.get();
}

} // namespace gbbs