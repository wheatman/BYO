

// Usage:
// numactl -i all ./run_unweighted -src 10012 -s -m -rounds 3 twitter_SJ
// flags:
//   required:
//     -src: the source to compute the BFS from
//   optional:
//     -rounds : the number of times to run the algorithm
//     -c : indicate that the graph is compressed
//     -m : indicate that the graph should be mmap'd
//     -s : indicate that the graph is symmetric
//     -d : dump the output arrays to files, useful for debugging

#include <functional>

#ifdef HOMEGROWN
#define PARLAY 1
#endif

#include "SSTGraph/SparseMatrix.hpp"
#include "gbbs/bridge.h"

template <template <class W> class vertex_type, class W>
struct symmetric_SSTGraph_graph {
  using vertex = vertex_type<W>;
  using weight_type = W;
  static constexpr bool binary = std::is_same_v<gbbs::empty, W>;
  using vertex_weight_type = double;
  using edge_type = typename vertex::edge_type;

  size_t N() const { return nodes.num_nodes(); }

  auto degree(size_t i) const { return nodes.getDegree(i); }

  template <class F> void map_neighbors(size_t i, F f) const {
    if constexpr (binary) {
      W empty_weight = W();
      nodes.map_neighbors_no_early_exit(
          i, [&](auto src, auto dest) { return f(src, dest, empty_weight); },
          nullptr, false);
    } else {
      nodes.map_neighbors_no_early_exit<0>(i, f, nullptr, false);
    }
  }

  template <class F> void parallel_map_neighbors(size_t i, F f) const {
    if constexpr (binary) {
      W empty_weight = W();
      nodes.map_neighbors_no_early_exit(
          i, [&](auto src, auto dest) { return f(src, dest, empty_weight); },
          nullptr, true);
    } else {
      nodes.map_neighbors_no_early_exit<0>(i, f, nullptr, true);
    }
  }

  template <class F> void map_neighbors_early_exit(size_t i, F f) const {
    if constexpr (binary) {
      W empty_weight = W();
      nodes.map_neighbors_early_exit(
          i, [&](auto src, auto dest) { return f(src, dest, empty_weight); },
          nullptr, false);
    } else {
      nodes.map_neighbors_early_exit<0>(i, f, nullptr, false);
    }
  }

  template <class F>
  void parallel_map_neighbors_early_exit(size_t i, F f) const {
    if constexpr (binary) {
      W empty_weight = W();
      nodes.map_neighbors_early_exit(
          i, [&](auto src, auto dest) { return f(src, dest, empty_weight); },
          nullptr, true);
    } else {
      nodes.map_neighbors_early_exit<0>(i, f, nullptr, true);
    }
  }

  void insert_sorted_batch(std::tuple<uint32_t, uint32_t> *es, size_t n) {
    nodes.insert_batch(es, n, true);
  }

  void remove_sorted_batch(std::tuple<uint32_t, uint32_t> *es, size_t n) {
    nodes.remove_batch(es, n, true);
  }

  // ======================= Constructors and fields  ========================
  symmetric_SSTGraph_graph() : vertex_weights(nullptr) {}

  // for now use the same constructor as the exsisting one for ease, but
  // probably write a more general constuctor for everyone to use later
  // TODO(wheatman) figure out what to do with _deletion_fn, probably should
  // just use unique pointer
  symmetric_SSTGraph_graph(auto *v_data, size_t n, size_t m,
                           std::function<void()> _deletion_fn, edge_type *_e0,
                           vertex_weight_type *_vertex_weights = nullptr)
      : nodes(n, n), vertex_weights(_vertex_weights),
        deletion_fn(_deletion_fn) {

    gbbs::parallel_for(0, n, [&](uint32_t i) {
      for (size_t j = v_data[i].offset; j < v_data[i].offset + v_data[i].degree;
           j++) {
        if constexpr (binary) {
          nodes.insert(i, std::get<0>(_e0[j]));
        } else {
          nodes.insert(i, std::get<0>(_e0[j]), std::get<1>(_e0[j]));
        }
      }
    });
  }

  // Move constructor
  symmetric_SSTGraph_graph(symmetric_SSTGraph_graph &&other) noexcept {
    nodes = other.nodes;
    vertex_weights = other.vertex_weights;
    other.vertex_weights = nullptr;
    deletion_fn = std::move(other.deletion_fn);
  }

  // Move assignment
  symmetric_SSTGraph_graph &
  operator=(symmetric_SSTGraph_graph &&other) noexcept {
    nodes = other.nodes;
    vertex_weights = other.vertex_weights;
    other.vertex_weights = nullptr;
    deletion_fn = std::move(other.deletion_fn);
  }

  // Copy constructor
  symmetric_SSTGraph_graph(const symmetric_SSTGraph_graph &other)
      : nodes(other.nodes) {
    vertex_weights = nullptr;
    if (other.vertex_weights != nullptr) {
      vertex_weights = gbbs::new_array_no_init<vertex_weight_type>(N());
      gbbs::parallel_for(0, N(), [&](size_t i) {
        vertex_weights[i] = other.vertex_weights[i];
      });
    }
    deletion_fn = [&]() {
      if (vertex_weights != nullptr) {
        gbbs::free_array(vertex_weights, N());
      }
    };
  }

  ~symmetric_SSTGraph_graph() { deletion_fn(); }

  // Graph Data
  std::conditional<binary, SSTGraph::SparseMatrixV<true>,
                   SSTGraph::SparseMatrixV<true, weight_type>>::type nodes;
  vertex_weight_type *vertex_weights;
  // called to delete the graph
  std::function<void()> deletion_fn;
};

#include "../run_unweighted.h"

using graph_impl =
    symmetric_SSTGraph_graph<gbbs::symmetric_vertex, gbbs::empty>;

using graph_api = gbbs::full_api;

using graph_t = gbbs::Graph<graph_impl, /* symmetric */ true, graph_api>;

int main(int argc, char *argv[]) {
  gbbs::commandLine P(argc, argv, " [-s] <inFile>");
  char *iFile = P.getArgument(0);
  bool symmetric = P.getOptionValue("-s");
  bool compressed = P.getOptionValue("-c");
  bool binary = P.getOptionValue("-b");
  bool mmap = P.getOptionValue("-m");
  gbbs::run_all_options options;
  options.dump = P.getOptionValue("-d");
  options.rounds = P.getOptionLongValue("-rounds", 3);
  options.src = static_cast<gbbs::uintE>(P.getOptionLongValue("-src", 0));

  std::cout << "### Graph: " << iFile << std::endl;
  if (compressed) {
    std::cerr << "is always compressed, but reads in uncompressed files\n";
    return -1;
  } else {
    if (symmetric) {
      auto G = gbbs::gbbs_io::read_unweighted_symmetric_graph<graph_t>(
          iFile, mmap, binary);
      run_all<true>(G, options);
    } else {
      std::cerr << "does not support directed graphs yet\n";
      return -1;
    }
  }
  return 1;
}