#pragma once

#include "gbbs/gbbs.h"
#include "reducer.hpp"

namespace gbbs {

template <class Graph> size_t Test_Map(const Graph &G) {
  Reducer_sum<size_t> count;

  parallel_for(0, G.N(), [&count, &G](size_t i) {
    size_t local_count = 0;
    G.map_neighbors(
        i, [&](auto src, auto dest, auto &val) { local_count += dest; });
    count += local_count;
  });
  return count.get();
}

} // namespace gbbs