// VeloGraphX adapter for the BYO graph-container benchmark framework.
//
// This adapter intentionally exposes VeloGraphX's mutable storage through
// BYO's graph-container API without using VeloGraphX's incremental analytics or
// adaptive execution policy. This keeps the comparison focused on the graph
// representation itself.

#include <cstddef>
#include <cstdint>
#include <functional>
#include <tuple>
#include <utility>
#include <vector>

#ifdef HOMEGROWN
#define PARLAY 1
#endif

#include "gbbs/bridge.h"
#include "velographx/storage/dynamic_graph.hpp"

#include "../run_unweighted.h"

template <class W>
struct VeloGraphX_graph {
  using vertex_weight_type = double;
  using weight_type = W;
  using edge_type = std::tuple<gbbs::uintE, W>;

  static constexpr bool insertable = true;
  static constexpr bool supports_get_memory_size = true;

  size_t N() const { return graph.vertex_count(); }
  size_t M() const { return graph.edge_count_directed(); }

  gbbs::uintE degree(size_t i) const {
    gbbs::uintE count = 0;
    graph.for_each_neighbor(static_cast<velographx::VertexId>(i),
                            [&](velographx::VertexId) { ++count; });
    return count;
  }

  template <class F>
  void map_neighbors(size_t i, F f) const {
    const auto src = static_cast<gbbs::uintE>(i);
    graph.for_each_neighbor(static_cast<velographx::VertexId>(i),
                            [&](velographx::VertexId dst) {
                              f(src, static_cast<gbbs::uintE>(dst), W{});
                            });
  }

  template <class F>
  void map_neighbors_early_exit(size_t i, F f) const {
    const auto src = static_cast<gbbs::uintE>(i);
    bool keep_going = true;
    graph.for_each_neighbor(static_cast<velographx::VertexId>(i),
                            [&](velographx::VertexId dst) {
                              if (keep_going &&
                                  f(src, static_cast<gbbs::uintE>(dst), W{})) {
                                keep_going = false;
                              }
                            });
  }

  void insert_sorted_batch(std::tuple<uint32_t, uint32_t>* es, size_t n) {
    velographx::UpdateBatch batch;
    batch.updates.reserve(n);
    for (size_t i = 0; i < n; ++i) {
      batch.add(static_cast<velographx::VertexId>(std::get<0>(es[i])),
                static_cast<velographx::VertexId>(std::get<1>(es[i])));
    }
    graph.apply(batch);
  }

  void remove_sorted_batch(std::tuple<uint32_t, uint32_t>* es, size_t n) {
    velographx::UpdateBatch batch;
    batch.updates.reserve(n);
    for (size_t i = 0; i < n; ++i) {
      batch.remove(static_cast<velographx::VertexId>(std::get<0>(es[i])),
                   static_cast<velographx::VertexId>(std::get<1>(es[i])));
    }
    graph.apply(batch);
  }

  VeloGraphX_graph() = default;

  VeloGraphX_graph(auto* v_data, size_t n, size_t m,
                   std::function<void()> cleanup, edge_type* e0,
                   vertex_weight_type* vertex_weights = nullptr)
      : graph(n, /* directed = */ true),
        deletion_fn(std::move(cleanup)),
        vertex_weights(vertex_weights) {
    std::vector<std::pair<velographx::VertexId, velographx::VertexId>> edges;
    edges.reserve(m);
    for (size_t u = 0; u < n; ++u) {
      const size_t begin = v_data[u].offset;
      const size_t end = begin + v_data[u].degree;
      for (size_t j = begin; j < end; ++j) {
        edges.emplace_back(static_cast<velographx::VertexId>(u),
                           static_cast<velographx::VertexId>(std::get<0>(e0[j])));
      }
    }
    graph.bulk_load_edges(edges);
  }

  VeloGraphX_graph(VeloGraphX_graph&& other) noexcept
      : graph(std::move(other.graph)),
        deletion_fn(std::move(other.deletion_fn)),
        vertex_weights(other.vertex_weights) {
    other.deletion_fn = []() {};
    other.vertex_weights = nullptr;
  }

  VeloGraphX_graph& operator=(VeloGraphX_graph&& other) noexcept {
    if (this != &other) {
      deletion_fn();
      graph = std::move(other.graph);
      deletion_fn = std::move(other.deletion_fn);
      vertex_weights = other.vertex_weights;
      other.deletion_fn = []() {};
      other.vertex_weights = nullptr;
    }
    return *this;
  }

  VeloGraphX_graph(const VeloGraphX_graph& other)
      : graph(other.graph), deletion_fn([]() {}), vertex_weights(nullptr) {}

  VeloGraphX_graph& operator=(const VeloGraphX_graph& other) {
    if (this != &other) {
      deletion_fn();
      graph = other.graph;
      deletion_fn = []() {};
      vertex_weights = nullptr;
    }
    return *this;
  }

  ~VeloGraphX_graph() { deletion_fn(); }

  size_t get_memory_size() { return graph.storage_bytes(); }

  velographx::DynamicGraph graph;
  std::function<void()> deletion_fn = []() {};
  vertex_weight_type* vertex_weights = nullptr;
};

using graph_impl = VeloGraphX_graph<gbbs::empty>;

// VeloGraphX currently exposes serial per-row iteration. BYO still parallelizes
// work across vertices; declaring no_parallel_map avoids claiming a native
// parallel neighbor-map primitive that the container does not provide.
using graph_api = gbbs::no_parallel_map;
using graph_t = gbbs::Graph<graph_impl, /* symmetric = */ true, graph_api>;

int main(int argc, char* argv[]) {
  gbbs::commandLine P(argc, argv, " [-s] <inFile>");
  char* iFile = P.getArgument(0);
  bool symmetric = P.getOptionValue("-s");
  bool compressed = P.getOptionValue("-c");
  bool binary = P.getOptionValue("-b");
  bool mmap = P.getOptionValue("-m");

  gbbs::run_all_options options;
  options.dump = P.getOptionValue("-d");
  options.rounds = P.getOptionLongValue("-rounds", 3);
  options.max_batch =
      static_cast<size_t>(P.getOptionLongValue("-max_batch", 1000000));
  options.src = static_cast<gbbs::uintE>(P.getOptionLongValue("-src", 0));
  options.inserts = P.getOptionValue("-i");

  std::cout << "### Graph: " << iFile << std::endl;
  if (compressed) {
    std::cerr << "VeloGraphX BYO adapter currently reads uncompressed inputs\n";
    return -1;
  }
  if (!symmetric) {
    std::cerr << "VeloGraphX BYO adapter currently benchmarks symmetric graphs only\n";
    return -1;
  }

  auto G = gbbs::gbbs_io::read_unweighted_symmetric_graph<graph_t>(
      iFile, mmap, binary);
  std::cout << "total bytes used = " << G.get_memory_size() << "\n";
  run_all(G, options);
  return 0;
}
